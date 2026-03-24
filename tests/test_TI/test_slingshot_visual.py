"""
Slingshot 功能测试脚本 - 使用模拟数据验证绘图功能

运行此脚本将：
1. 生成模拟精子发生轨迹数据
2. 运行 Slingshot 轨迹推断
3. 生成并保存交互式可视化 HTML 文件
"""
import numpy as np
from anndata import AnnData


def generate_test_data_spermatogenesis(n_cells=520, random_state=42):
    """
    生成模拟精子发生轨迹的测试数据

    模拟细胞分化路径：
    SSC → Undiff.SPG → Diff.SPG → preL → L → Z → preP → P → D → SPC → S1 → S2 → S3 → S4
    """
    np.random.seed(random_state)

    cell_types = [
        "SSC", "Undiff.SPG", "Diff.SPG", "preL",
        "L", "Z", "preP", "P", "D", "SPC",
        "S1", "S2", "S3", "S4"
    ]

    n_cell_types = len(cell_types)
    cells_per_type = [n_cells // n_cell_types] * n_cell_types
    cells_per_type[0] += n_cells - sum(cells_per_type)

    all_cells = []
    all_labels = []

    for i, (ct, n) in enumerate(zip(cell_types, cells_per_type)):
        t = i / n_cell_types * np.pi
        base_x = np.cos(t) * 5
        base_y = np.sin(t) * 5

        cells_x = base_x + np.random.randn(n) * 0.8
        cells_y = base_y + np.random.randn(n) * 0.8

        all_cells.append(np.column_stack([cells_x, cells_y]))
        all_labels.extend([ct] * n)

    X_umap = np.vstack(all_cells)

    adata = AnnData(
        X=np.random.randn(n_cells, 100),
        obs={"celltype": all_labels},
        obsm={"X_umap": X_umap}
    )

    predefined_tree = {
        "SSC": ["Undiff.SPG"],
        "Undiff.SPG": ["Diff.SPG"],
        "Diff.SPG": ["preL"],
        "preL": ["L"],
        "L": ["Z"],
        "Z": ["preP"],
        "preP": ["P"],
        "P": ["D"],
        "D": ["SPC"],
        "SPC": ["S1"],
        "S1": ["S2"],
        "S2": ["S3"],
        "S3": ["S4"]
    }

    adata.uns["predefined_tree"] = predefined_tree
    return adata


def generate_test_data_bifurcation(n_cells=600, random_state=42):
    """
    生成具有分支的测试数据
    """
    np.random.seed(random_state)

    cell_types = [
        "SSC", "Undiff.SPG", "Diff.SPG", "preL",
        "L", "Z", "L2", "L3", "preP", "P", "D", "SPC"
    ]

    n_cell_types = len(cell_types)
    cells_per_type = [n_cells // n_cell_types] * n_cell_types
    cells_per_type[0] += n_cells - sum(cells_per_type)

    all_cells = []
    all_labels = []

    for i, (ct, n) in enumerate(zip(cell_types, cells_per_type)):
        if ct in ["L", "Z"]:
            base_x, base_y = 0, 0
        elif ct in ["L2", "L3", "preP", "P", "D", "SPC"]:
            idx = cell_types.index(ct) - cell_types.index("L2")
            t = idx / 6 * np.pi / 3
            base_x = np.cos(t) * 4 + 2
            base_y = np.sin(t) * 4 + 2
        else:
            idx = cell_types.index(ct)
            t = idx / 4 * np.pi / 2
            base_x = -np.cos(t) * 5
            base_y = -np.sin(t) * 5

        cells_x = base_x + np.random.randn(n) * 0.6
        cells_y = base_y + np.random.randn(n) * 0.6

        all_cells.append(np.column_stack([cells_x, cells_y]))
        all_labels.extend([ct] * n)

    X_umap = np.vstack(all_cells)

    adata = AnnData(
        X=np.random.randn(n_cells, 100),
        obs={"celltype": all_labels},
        obsm={"X_umap": X_umap}
    )

    predefined_tree = {
        "SSC": ["Undiff.SPG"],
        "Undiff.SPG": ["Diff.SPG"],
        "Diff.SPG": ["preL"],
        "preL": ["L"],
        "L": ["Z", "L2"],
        "Z": ["preP"],
        "L2": ["L3"],
        "L3": ["preP"],
        "preP": ["P"],
        "P": ["D"],
        "D": ["SPC"]
    }

    adata.uns["predefined_tree"] = predefined_tree
    return adata


def run_visual_test(test_name="spermatogenesis", output_dir="tests/test_TI/output"):
    """
    运行可视化测试

    Parameters
    ----------
    test_name : str
        测试名称 ("spermatogenesis" 或 "bifurcation")
    output_dir : str
        输出目录
    """
    import os
    from hic_basic.TI.slingshot import Slingshot
    from hic_basic.TI.slingshot.plot import plot_slingshot, plot_clusters, plot_curves

    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Slingshot Visual Test: {test_name}")
    print(f"{'='*60}\n")

    # 1. 生成测试数据
    if test_name == "spermatogenesis":
        adata = generate_test_data_spermatogenesis(n_cells=520)
    elif test_name == "bifurcation":
        adata = generate_test_data_bifurcation(n_cells=600)
    else:
        raise ValueError(f"Unknown test_name: {test_name}")

    print(f"Generated {adata.n_obs} cells")
    print(f"Cell types: {adata.obs['celltype'].unique().tolist()}")

    # 2. 初始化 Slingshot
    predefined_tree = adata.uns["predefined_tree"]

    slingshot_obj = Slingshot(
        data=adata,
        celltype_key="celltype",
        obsm_key="X_umap",
        predefined_tree=predefined_tree,
        start_node="SSC",
        debug_level=None
    )

    print(f"\nSlingshot initialized")
    print(f"Number of clusters: {slingshot_obj.num_clusters}")
    print(f"Start node: {slingshot_obj.start_node}")

    # 3. 拟合
    print("\nFitting Slingshot...")
    slingshot_obj.fit(num_epochs=10)
    print(f"Fit complete. Number of lineages: {len(slingshot_obj.lineages)}")

    # 4. 打印谱系信息
    print("\nLineages:")
    slingshot_obj.list_lineages()

    # 5. 生成可视化
    print("\nGenerating visualizations...")

    # 主图：pseudotime + curves
    output_path_main = f"{output_dir}/slingshot_{test_name}_main.png"
    plot_slingshot(
        slingshot_obj,
        title=f"Slingshot Trajectory - {test_name.title()}",
        save_path=output_path_main,
        show=False
    )
    print(f"  ✓ Main plot saved: {output_path_main}")

    # 聚类图
    output_path_clusters = f"{output_dir}/slingshot_{test_name}_clusters.png"
    plot_clusters(
        slingshot_obj,
        title=f"Cell Clusters - {test_name.title()}",
        save_path=output_path_clusters,
        show=False
    )
    print(f"  ✓ Clusters plot saved: {output_path_clusters}")

    # 曲线图
    output_path_curves = f"{output_dir}/slingshot_{test_name}_curves.png"
    plot_curves(
        slingshot_obj,
        title=f"Principal Curves - {test_name.title()}",
        save_path=output_path_curves,
        show=False
    )
    print(f"  ✓ Curves plot saved: {output_path_curves}")

    # 6. 验证 pseudotime
    pseudotime = slingshot_obj.unified_pseudotime
    print(f"\nPseudotime statistics:")
    print(f"  Min: {pseudotime.min():.4f}")
    print(f"  Max: {pseudotime.max():.4f}")
    print(f"  Mean: {pseudotime.mean():.4f}")
    print(f"  Std: {pseudotime.std():.4f}")

    print(f"\n{'='*60}")
    print("Visual test complete!")
    print(f"Output directory: {output_dir}/")
    print(f"{'='*60}\n")

    return slingshot_obj


if __name__ == "__main__":
    import sys

    test_name = sys.argv[1] if len(sys.argv) > 1 else "spermatogenesis"
    run_visual_test(test_name)
