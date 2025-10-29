#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, linregress
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import argparse
import warnings
import os 


def safe_pearsonr(x, y):

    x = np.asarray(x)
    y = np.asarray(y)
    if x.size == 0 or y.size == 0 or np.all(x == x.flat[0]) or np.all(y == y.flat[0]):
        return np.nan, np.nan 

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            pearson_constant_warning = getattr(pearsonr, 'ConstantInputWarning', RuntimeWarning)
            r, p = pearsonr(x, y)
            if any(isinstance(wrn.message, pearson_constant_warning) or \
                   "constant" in str(wrn.message).lower() for wrn in w):
                return np.nan, np.nan
        except ValueError:
             return np.nan, np.nan
    return r, p

def compute_t_proj(df, x0, y0, x1, y1):
    df = df.copy()
    dx, dy = x1 - x0, y1 - y0
    d_len_sq = dx**2 + dy**2
    if d_len_sq == 0:
        df['t_proj'] = np.nan
        return df
    d_len = np.sqrt(d_len_sq)
    # Ensure 'x' and 'y' columns are numeric
    df['x'] = pd.to_numeric(df['x'], errors='coerce')
    df['y'] = pd.to_numeric(df['y'], errors='coerce')
    df['t_proj'] = ((df['x'] - x0) * dx + (df['y'] - y0) * dy) / d_len
    return df

def transform_coordinates(df_in, x0, y0, x1, y1, flip_axis='none'):
    df = df_in.copy()

    if not all(col in df.columns for col in ['x', 'y']):
        print("警告: transform_coordinates 需要 DataFrame 包含 'x' 和 'y' 列。")
        df['x_transformed'] = np.nan
        df['y_transformed'] = np.nan
        return df, (0,0), (0,0)

    dx = x1 - x0
    dy = y1 - y0
    length_p0p1 = np.hypot(dx, dy)


    df['x_translated'] = df['x'] - x0
    df['y_translated'] = df['y'] - y0

    p0_transformed_final = (0.0, 0.0) 
    p1_x_final = length_p0p1
    p1_y_final = 0.0

    if length_p0p1 == 0: 
        df['x_transformed'] = df['x_translated']
        df['y_transformed'] = df['y_translated']
    else:
        cos_theta = dx / length_p0p1
        sin_theta = dy / length_p0p1

        df['x_transformed'] = df['x_translated'] * cos_theta + df['y_translated'] * sin_theta
        df['y_transformed'] = -df['x_translated'] * sin_theta + df['y_translated'] * cos_theta

    if 'x' in flip_axis:
        df['x_transformed'] = -df['x_transformed']
        p1_x_final = -p1_x_final 
        if p0_transformed_final[0] != 0: 
            p0_transformed_final = (-p0_transformed_final[0], p0_transformed_final[1])


    if 'y' in flip_axis:
        df['y_transformed'] = -df['y_transformed']
        if p0_transformed_final[1] != 0: 
             p0_transformed_final = (p0_transformed_final[0], -p0_transformed_final[1])


    p1_transformed_final = (p1_x_final, p1_y_final)

    df.drop(columns=['x_translated', 'y_translated'], inplace=True, errors='ignore')
    return df, p0_transformed_final, p1_transformed_final


def plot_spatial_expression(gem_df_target_gene, gem_df_other_genes_unique_coords,
                            x0, y0, x1, y1, gene_id, output_base_path):
    plt.figure(figsize=(10, 8))

    if gem_df_other_genes_unique_coords is not None and not gem_df_other_genes_unique_coords.empty:
        plt.scatter(gem_df_other_genes_unique_coords['x'], gem_df_other_genes_unique_coords['y'],
                    s=5, color='lightgray', alpha=0.5, label='Other genes (unique positions)')

    target_gene_plot_df = gem_df_target_gene.copy() 
    if not target_gene_plot_df.empty:
        target_gene_plot_df['MIDCount'] = pd.to_numeric(target_gene_plot_df['MIDCount'], errors='coerce')
        target_gene_plot_df.dropna(subset=['MIDCount', 'x', 'y'], inplace=True)

        if not target_gene_plot_df.empty:
            min_val = target_gene_plot_df['MIDCount'].min()
            max_val = target_gene_plot_df['MIDCount'].max()
            norm = None
            if min_val > 0 and max_val > min_val :
                norm = mcolors.LogNorm(vmin=min_val, vmax=max_val)
            elif max_val > min_val:
                norm = mcolors.Normalize(vmin=min_val, vmax=max_val)
            elif min_val == max_val :
                norm = mcolors.Normalize(vmin=min_val - 0.5 if min_val > 0 else min_val,
                                         vmax=max_val + 0.5 if max_val > 0 else max_val +1)
                if norm.vmin >= norm.vmax:
                    norm.vmin = 0
                    norm.vmax = 1 if max_val == 0 else max_val * 2 if max_val > 0 else 1

            if norm is not None:
                scatter = plt.scatter(target_gene_plot_df['x'], target_gene_plot_df['y'],
                                      s=5, c=target_gene_plot_df['MIDCount'], cmap='viridis', 
                                      norm=norm, label=f'Gene {gene_id} expression')
                cbar = plt.colorbar(scatter, label='MIDCount (Expression Level)')
            else: 
                 plt.scatter(target_gene_plot_df['x'], target_gene_plot_df['y'],
                                      s=5, color='blue', label=f'Gene {gene_id} expression (uniform/no scale)')
        else:
            print(f"警告: 基因 {gene_id} 没有有效的数值型 MIDCount 用于常规空间绘图。")


    plt.plot([x0, x1], [y0, y1], color='red', linestyle='--', linewidth=2, label='Projection Line (P0-P1)')
    plt.scatter([x0, x1], [y0, y1], color='red', s=50, marker='o')

    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title(f'Spatial Expression of Gene: {gene_id}')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.axis('equal')

    output_path_png = output_base_path + ".png"
    output_path_pdf = output_base_path + ".pdf"

    try:
        plt.savefig(output_path_png, dpi=300, bbox_inches='tight')
        print(f" {output_path_png}")
        plt.savefig(output_path_pdf, dpi=300, bbox_inches='tight')
        print(f"{output_path_pdf}")
    except Exception as e:
        print(f"{e}")
    plt.close()

def plot_transformed_spatial_expression(gem_df_target_gene, gem_df_other_genes_unique_coords,
                                        x0, y0, x1, y1, gene_id, flip_axis_param, output_base_path_transform):
    target_gene_transformed_df, p0_new, p1_new = transform_coordinates(gem_df_target_gene.copy(), x0, y0, x1, y1, flip_axis=flip_axis_param)
    other_genes_transformed_df = None
    if gem_df_other_genes_unique_coords is not None and not gem_df_other_genes_unique_coords.empty:
        other_genes_transformed_df, _, _ = transform_coordinates(gem_df_other_genes_unique_coords.copy(), x0, y0, x1, y1, flip_axis=flip_axis_param)

    plt.figure(figsize=(10, 8))
    if other_genes_transformed_df is not None and not other_genes_transformed_df.empty:
        plt.scatter(other_genes_transformed_df['x_transformed'], other_genes_transformed_df['y_transformed'],
                    s=5, color='lightgray', alpha=0.5, label='Other genes (transformed)')
    if not target_gene_transformed_df.empty:
        target_gene_transformed_df['MIDCount'] = pd.to_numeric(target_gene_transformed_df['MIDCount'], errors='coerce')
        plot_df = target_gene_transformed_df.dropna(subset=['x_transformed', 'y_transformed', 'MIDCount'])

        if not plot_df.empty:
            min_val = plot_df['MIDCount'].min()
            max_val = plot_df['MIDCount'].max()
            norm = None
            if min_val > 0 and max_val > min_val :
                norm = mcolors.LogNorm(vmin=min_val, vmax=max_val)
            elif max_val > min_val:
                norm = mcolors.Normalize(vmin=min_val, vmax=max_val)
            elif min_val == max_val :
                norm = mcolors.Normalize(vmin=min_val - 0.5 if min_val > 0 else min_val,
                                         vmax=max_val + 0.5 if max_val > 0 else max_val +1)
                if norm.vmin >= norm.vmax:
                    norm.vmin = 0
                    norm.vmax = 1 if max_val == 0 else max_val * 2 if max_val > 0 else 1

            if norm is not None:
                scatter = plt.scatter(plot_df['x_transformed'], plot_df['y_transformed'],
                                      s=5, c=plot_df['MIDCount'], cmap='viridis', # Changed cmap to coolwarm
                                      norm=norm, label=f'Gene {gene_id} expression (transformed)')
                cbar = plt.colorbar(scatter, label='MIDCount (Expression Level)')
            else:
                 plt.scatter(plot_df['x_transformed'], plot_df['y_transformed'],
                                      s=5, color='blue', label=f'Gene {gene_id} expression (transformed, uniform/no scale)')
        else:
            print("Error")


    plt.plot([p0_new[0], p1_new[0]], [p0_new[1], p1_new[1]],
             color='red', linestyle='--', linewidth=2, label='Projection Line (transformed P0-P1)')
    plt.scatter([p0_new[0], p1_new[0]], [p0_new[1], p1_new[1]], color='red', s=50, marker='o')


    plt.xlabel('X Coordinate (Transformed)')
    plt.ylabel('Y Coordinate (Transformed)')
    plt.title(f'Transformed Spatial Expression of Gene: {gene_id} (P0-P1 aligned to X-axis, Flip: {flip_axis_param})')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.axis('equal')

    output_path_png = output_base_path_transform + ".png"
    output_path_pdf = output_base_path_transform + ".pdf"

    try:
        plt.savefig(output_path_png, dpi=300, bbox_inches='tight')
        print(f"{output_path_png}")
        plt.savefig(output_path_pdf, dpi=300, bbox_inches='tight')
        print(f"{output_path_pdf}")
    except Exception as e:
        print(f"{e}")
    plt.close()


def plot_expression_projection(gem_df_target_gene, gene_id, min_points, output_base_path):
    plt.figure(figsize=(10, 6))

    plot_target_df = gem_df_target_gene.copy()
    plot_target_df['t_proj'] = pd.to_numeric(plot_target_df['t_proj'], errors='coerce')
    plot_target_df['MIDCount'] = pd.to_numeric(plot_target_df['MIDCount'], errors='coerce')
    plot_data = plot_target_df[['t_proj', 'MIDCount']].dropna()

    num_valid_points = len(plot_data)

    plt.scatter(plot_data['t_proj'], plot_data['MIDCount'], alpha=0.6, label='Data points')

    if plot_data.empty or num_valid_points < min_points:
        print(f"Error: (N={num_valid_points}) < ({min_points})")
    else:
        try:
            slope, intercept, r_value, p_value, std_err = linregress(plot_data['t_proj'], plot_data['MIDCount'])

            x_plot_min, x_plot_max = plot_data['t_proj'].min(), plot_data['t_proj'].max()
            if x_plot_min == x_plot_max:
                x_vals_reg = np.array([x_plot_min -1, x_plot_max +1])
            else:
                x_vals_reg = np.array([x_plot_min, x_plot_max])

            y_vals_reg = intercept + slope * x_vals_reg
            plt.plot(x_vals_reg, y_vals_reg, color='red', label=f'Linear Regression\ny = {slope:.2f}x + {intercept:.2f}')

            text_str = f'N = {num_valid_points}\nPearson R: {r_value:.3f}\nP-value: {p_value:.3g}'
            plt.text(0.05, 0.95, text_str, transform=plt.gca().transAxes, fontsize=10,
                     verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))
        except ValueError as e:
            error_text = f' (N={num_valid_points})'
            print(f"Error: {e}. ")
            plt.text(0.5, 0.5, error_text,
                     ha='center', va='center', transform=plt.gca().transAxes,
                     fontsize=12, color='red')

    plt.xlabel('Projection on Line (t_proj)')
    plt.ylabel('MIDCount (Expression Level)')
    plt.title(f'Expression vs. Projection for Gene: {gene_id}')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.7)

    output_path_png = output_base_path + ".png"
    output_path_pdf = output_base_path + ".pdf"
    try:
        plt.savefig(output_path_png, dpi=300, bbox_inches='tight')
        print(f"Error: {output_path_png}")
        plt.savefig(output_path_pdf, dpi=300, bbox_inches='tight')
        print(f"Error: {output_path_pdf}")
    except Exception as e:
        print(f"Error: {e}")
    plt.close()

def parse_args():
    parser = argparse.ArgumentParser(
        description="为指定基因绘制空间表达图和沿直线方向的表达趋势图。"
    )
    parser.add_argument("gem_path", help="输入 GEM 文件路径 (TSV 格式，含 header)")
    parser.add_argument("--gene_to_plot", "-g", type=str, required=True, help="要绘制的目标基因的 geneID")
    parser.add_argument("--x0",   type=float, required=True, help="直线端点 P0 的 x 坐标")
    parser.add_argument("--y0",   type=float, required=True, help="直线端点 P0 的 y 坐标")
    parser.add_argument("--x1",   type=float, required=True, help="直线端点 P1 的 x 坐标")
    parser.add_argument("--y1",   type=float, required=True, help="直线端点 P1 的 y 坐标")
    parser.add_argument("--output_prefix", "-o", default="gene_plot",
                        help="输出图片文件的前缀 (默认: gene_plot)")
    parser.add_argument("--min-points", "-m", type=int, default=10,
                        help="进行回归分析和显示统计值所需的最少点数 (默认: 10)")
    parser.add_argument("--xf", type=float, help="过滤参考点的 x 坐标 (可选, 若使用过滤则为必需)")
    parser.add_argument("--yf", type=float, help="过滤参考点的 y 坐标 (可选, 若使用过滤则为必需)")
    parser.add_argument("--filter_action", choices=['remove', 'keep'],
                        help="过滤操作：'remove' 移除匹配点, 'keep' 保留匹配点 (可选, 若使用过滤则为必需)")
    parser.add_argument("--x_comparison", choices=['>', '<', '='],
                        help="x 坐标的比较操作符 (可选, 若使用过滤则为必需)")
    parser.add_argument("--y_comparison", choices=['>', '<', '='],
                        help="y 坐标的比较操作符 (可选, 若使用过滤则为必需)")

    # New argument for flipping transformed axes
    parser.add_argument("--transform_flip", choices=['x', 'y', 'xy', 'none'], default='none',
                        help="在坐标变换后翻转轴：'x', 'y', 'xy', 或 'none' (默认: none)")
    return parser.parse_args()

def main():
    args = parse_args()

    try:
        gem_full = pd.read_csv(args.gem_path, sep='\t', header=0, comment='#',
                               usecols=['geneID', 'x', 'y', 'MIDCount'])
    except FileNotFoundError:
        print(f"Error")
        return
    except ValueError as e:
        print(f"Error")
        return

    if gem_full.empty:
        print("Error")
        return
    print("Error")

    gem_full['x'] = pd.to_numeric(gem_full['x'], errors='coerce')
    gem_full['y'] = pd.to_numeric(gem_full['y'], errors='coerce')
    gem_full['MIDCount'] = pd.to_numeric(gem_full['MIDCount'], errors='coerce')
    gem_full.dropna(subset=['x', 'y', 'MIDCount'], inplace=True)
    if gem_full.empty:
        print("Error")
        return

    gem_filtered = gem_full.copy()
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
        initial_rows = len(gem_filtered)

        condition_x = None
        if args.x_comparison == '>': condition_x = (gem_filtered['x'] > args.xf)
        elif args.x_comparison == '<': condition_x = (gem_filtered['x'] < args.xf)
        elif args.x_comparison == '=': condition_x = (gem_filtered['x'] == args.xf)

        condition_y = None
        if args.y_comparison == '>': condition_y = (gem_filtered['y'] > args.yf)
        elif args.y_comparison == '<': condition_y = (gem_filtered['y'] < args.yf)
        elif args.y_comparison == '=': condition_y = (gem_filtered['y'] == args.yf)

    elif some_filtering_args_provided:
        print("Error")

    gem_projected = compute_t_proj(gem_filtered.copy(), args.x0, args.y0, args.x1, args.y1)
    gem_projected.dropna(subset=['t_proj'], inplace=True)
    if gem_projected.empty:
        print("Error")
        return

    target_gene_data = gem_projected[gem_projected['geneID'] == args.gene_to_plot].copy()
    other_genes_data = gem_projected[gem_projected['geneID'] != args.gene_to_plot].copy()

    if target_gene_data.empty:
        print(f"Error")
        return


    other_genes_unique_coords = None
    if not other_genes_data.empty:
        other_genes_unique_coords = other_genes_data[['x', 'y']].drop_duplicates().reset_index(drop=True)

    plot_spatial_expression(target_gene_data.copy(),
                            other_genes_unique_coords.copy() if other_genes_unique_coords is not None else None,
                            args.x0, args.y0, args.x1, args.y1,
                            args.gene_to_plot, f"{args.output_prefix}_spatial")

    plot_transformed_spatial_expression(target_gene_data.copy(),
                                        other_genes_unique_coords.copy() if other_genes_unique_coords is not None else None,
                                        args.x0, args.y0, args.x1, args.y1,
                                        args.gene_to_plot,
                                        args.transform_flip, 
                                        f"{args.output_prefix}_spatial_transform")

    plot_expression_projection(target_gene_data.copy(), args.gene_to_plot,
                               args.min_points,
                               f"{args.output_prefix}_projection")

    print(f"Done")

if __name__ == '__main__':
    main()