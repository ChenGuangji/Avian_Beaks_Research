#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import warnings

def safe_pearsonr(x, y):

    x = np.asarray(x)
    y = np.asarray(y)

    if x.size == 0 or y.size == 0 or np.all(x == x.flat[0]) or np.all(y == y.flat[0]):
        return np.nan, np.nan, 0 

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:

            pearson_constant_warning = getattr(pearsonr, 'ConstantInputWarning', RuntimeWarning) # Fallback for older scipy

            r, p = pearsonr(x, y)

            if any(isinstance(wrn.message, pearson_constant_warning) or \
                   "constant" in str(wrn.message).lower() for wrn in w):
                return np.nan, np.nan, x.size 
        except ValueError:
                 return np.nan, np.nan, 0
    return r, p, x.size 

def safe_spearmanr(x, y): 
    x = np.asarray(x)
    y = np.asarray(y)
    if x.size == 0 or y.size == 0 or np.all(x == x.flat[0]) or np.all(y == y.flat[0]):
        return np.nan, np.nan, 0

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always") 
        try:
            r, p = spearmanr(x, y)
            if any("constant" in str(wrn.message).lower() for wrn in w): 
                return np.nan, np.nan, x.size
        except ValueError: 
            return np.nan, np.nan, 0
    return r, p, x.size

def compute_t_proj(df, x0, y0, x1, y1):

    dx, dy = x1 - x0, y1 - y0
    d_len_sq = dx**2 + dy**2 
    if d_len_sq == 0: 
        df = df.copy()
        df['t_proj'] = np.nan
        return df
    d_len = np.sqrt(d_len_sq)
    df = df.copy()

    df['x'] = pd.to_numeric(df['x'], errors='coerce')
    df['y'] = pd.to_numeric(df['y'], errors='coerce')
    df['t_proj'] = ((df['x'] - x0) * dx + (df['y'] - y0) * dy) / d_len
    return df

def gene_corr(args_tuple):

    gene, subdf, min_points = args_tuple 


    if not isinstance(subdf, pd.DataFrame) or 't_proj' not in subdf.columns or 'MIDCount' not in subdf.columns:
        return gene, np.nan, np.nan, np.nan, 0 


    if len(subdf) < min_points: 
        return gene, np.nan, np.nan, np.nan, 0 

    t_proj_values = pd.to_numeric(subdf['t_proj'], errors='coerce')
    midcount_values = pd.to_numeric(subdf['MIDCount'], errors='coerce')

    
    temp_df = pd.DataFrame({'t_proj': t_proj_values, 'MIDCount': midcount_values}).dropna()

    num_valid_points_after_nan_removal = len(temp_df)

    if num_valid_points_after_nan_removal < min_points:
        return gene, np.nan, np.nan, np.nan, num_valid_points_after_nan_removal

    # Use the cleaned data for correlation
    corr, pval, _ = safe_pearsonr(temp_df['t_proj'], temp_df['MIDCount'])

    if np.isnan(corr):
        r_squared = np.nan
        return gene, np.nan, np.nan, r_squared, num_valid_points_after_nan_removal
    else:
        r_squared = corr**2

    return gene, corr, pval, r_squared, num_valid_points_after_nan_removal

def parse_args():
    parser = argparse.ArgumentParser(
        description="Pearson 相关"
    )
    parser.add_argument("gem_path", help="输入 GEM 文件路径 (TSV 格式，含 header)")
    parser.add_argument("--x0",   type=float, required=True, help="直线端点 P0 的 x 坐标")
    parser.add_argument("--y0",   type=float, required=True, help="直线端点 P0 的 y 坐标")
    parser.add_argument("--x1",   type=float, required=True, help="直线端点 P1 的 x 坐标")
    parser.add_argument("--y1",   type=float, required=True, help="直线端点 P1 的 y 坐标")
    parser.add_argument("--min-points", "-m", type=int, default=5,
                        help="每个基因至少包含的有效数值点对数 (默认: 5)")
    parser.add_argument("--Cores",    "-c", type=int, default=4,
                        help="并行进程数 (默认: 4)")
    parser.add_argument("--output",     "-o", default="gene_gradient_correlation.tsv",
                        help="输出 TSV 文件路径 (默认: gene_gradient_correlation.tsv)")
    parser.add_argument("--xf", type=float, help="过滤参考点的 x 坐标 (可选, 若使用过滤则为必需)")
    parser.add_argument("--yf", type=float, help="过滤参考点的 y 坐标 (可选, 若使用过滤则为必需)")
    parser.add_argument("--filter_action", choices=['remove', 'keep'],
                        help="过滤操作：'remove' 移除匹配点, 'keep' 保留匹配点 (可选, 若使用过滤则为必需)")
    parser.add_argument("--x_comparison", choices=['>', '<', '='],
                        help="x 坐标的比较操作符 (可选, 若使用过滤则为必需)")
    parser.add_argument("--y_comparison", choices=['>', '<', '='],
                        help="y 坐标的比较操作符 (可选, 若使用过滤则为必需)")
    return parser.parse_args()

def main():
    args = parse_args()

    try:
        gem = pd.read_csv(args.gem_path, sep='\t', header=0,comment='#',
                          usecols=['geneID','x','y','MIDCount'])
    except FileNotFoundError:
        print(f"GEM  {args.gem_path}")
        return
    except ValueError as e:
        print(f"('geneID','x','y','MIDCount') {e}")
        return

    if gem.empty:
        print("GEM 文件为空")
        return

    filtering_args_provided = (args.xf is not None and
                               args.yf is not None and
                               args.filter_action is not None and
                               args.x_comparison is not None and
                               args.y_comparison is not None)

    some_filtering_args_provided = (args.xf is not None or
                                    args.yf is not None or
                                    args.filter_action is not None or
                                    args.x_comparison is not None or
                                    args.y_comparison is not None)

    if filtering_args_provided:
        initial_rows = len(gem)

        gem['x'] = pd.to_numeric(gem['x'], errors='coerce')
        gem['y'] = pd.to_numeric(gem['y'], errors='coerce')
        gem.dropna(subset=['x', 'y'], inplace=True)

        condition_x = None
        if args.x_comparison == '>': condition_x = (gem['x'] > args.xf)
        elif args.x_comparison == '<': condition_x = (gem['x'] < args.xf)
        elif args.x_comparison == '=': condition_x = (gem['x'] == args.xf)

        condition_y = None
        if args.y_comparison == '>': condition_y = (gem['y'] > args.yf)
        elif args.y_comparison == '<': condition_y = (gem['y'] < args.yf)
        elif args.y_comparison == '=': condition_y = (gem['y'] == args.yf)

    elif some_filtering_args_provided:
        print("--xf, --yf, --filter_action, --x_comparison,  --y_comparison ")


    gem = compute_t_proj(gem, args.x0, args.y0, args.x1, args.y1)

    gem.dropna(subset=['t_proj'], inplace=True)
    if gem.empty:
        print("Error")
        return

    tasks = []
    for gene, df_group in gem.groupby('geneID'):
        if isinstance(df_group, pd.DataFrame) and not df_group.empty:
            tasks.append((gene, df_group.copy(), args.min_points))
        else:
            print(f"Error")


    if not tasks:
        print("Error")
        return

    results = []
    num_cores = max(1, args.Cores)
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        future_to_gene = {executor.submit(gene_corr, task_args): task_args[0] for task_args in tasks}
        for i, fut in enumerate(as_completed(future_to_gene)):
            gene_name = future_to_gene[fut]
            try:
                gene, corr, pval, r_squared, num_points = fut.result() 
                results.append((gene, corr, pval, r_squared, num_points)) 
            except Exception as exc:
                print(f"Error")
                results.append((gene_name, np.nan, np.nan, np.nan, 0)) 
            finally:
                print(f"Done: {i+1}/{len(tasks)} ({((i+1)/len(tasks)*100):.2f}%)", end="\r")

    if not results:
        print("Error")
        return

    res_df = pd.DataFrame(results, columns=['geneID','corr_coeff','p_value', 'R_squared', 'num_valid_points'])
    res_df['abs_corr_coeff'] = res_df['corr_coeff'].abs()

    res_df = res_df.sort_values(by=['p_value', 'abs_corr_coeff'], ascending=[True, False]).drop(columns=['abs_corr_coeff'])
    res_df = res_df.reset_index(drop=True)

    try:
        res_df.to_csv(args.output, sep='\t', index=False, na_rep='NaN')
    except Exception as e:
        print(f"Error")


if __name__ == '__main__':
    main()