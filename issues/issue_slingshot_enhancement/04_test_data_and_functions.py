"""
测试数据和绘图函数 - 用于开发后检查 Slingshot 模块

基于 spg 项目的真实数据生成的测试用例
"""

import numpy as np
from anndata import AnnData


def generate_test_data_spermatogenesis(n_cells=500, random_state=42):
    """
    生成模拟精子发生轨迹的测试数据
    
    模拟精子发生的细胞分化路径：
    SSC → Undiff.SPG → Diff.SPG → preL → L → Z → preP → P → D → SPC → S1 → S2 → S3 → S4
    
    参数:
        n_cells: 细胞数量
        random_state: 随机种子
    
    返回:
        adata: AnnData 对象，包含：
            - obs['celltype']: 细胞类型标签（字符串）
            - obsm['X_umap']: 2D UMAP 嵌入
            - uns['predefined_tree']: 预定义的树结构
    """
    np.random.seed(random_state)
    
    # 定义精子发生细胞类型（简化版）
    cell_types = [
        "SSC",           # 0: 精原干细胞
        "Undiff.SPG",    # 1: 未分化精母细胞
        "Diff.SPG",      # 2: 分化精母细胞
        "preL",          # 3: 前细线期
        "L",             # 4: 细线期
        "Z",             # 5: 偶线期
        "preP",          # 6: 前粗线期
        "P",             # 7: 粗线期
        "D",             # 8: 双线期
        "SPC",           # 9: 精细胞
        "S1",            # 10: 精子 1
        "S2",            # 11: 精子 2
        "S3",            # 12: 精子 3
        "S4"             # 13: 精子 4
    ]
    
    n_cell_types = len(cell_types)
    cells_per_type = n_cells // n_cell_types
    
    # 生成 UMAP 嵌入（沿轨迹的弧形结构）
    all_cells = []
    all_labels = []
    
    for i, ct in enumerate(cell_types):
        # 沿弧形的参数位置
        t = i / n_cell_types * np.pi  # 0 到 π
        base_x = np.cos(t) * 5
        base_y = np.sin(t) * 5
        
        # 添加细胞类型的变异
        cells_x = base_x + np.random.randn(cells_per_type) * 0.8
        cells_y = base_y + np.random.randn(cells_per_type) * 0.8
        
        all_cells.append(np.column_stack([cells_x, cells_y]))
        all_labels.extend([ct] * cells_per_type)
    
    X_umap = np.vstack(all_cells)
    
    # 创建 AnnData 对象
    adata = AnnData(
        X=np.random.randn(n_cells, 100),  # 假基因表达数据
        obs={'celltype': all_labels},
        obsm={'X_umap': X_umap}
    )
    
    # 定义预定义的树结构
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
    
    adata.uns['predefined_tree'] = predefined_tree
    
    return adata


def generate_test_data_bifurcation(n_cells=600, random_state=42):
    """
    生成具有分支的测试数据（如 L/Z 分支）
    
    模拟存在分支的分化路径：
    preL → L → L2 → L3 → preP
           ↘ Z ↗
    
    参数:
        n_cells: 细胞数量
        random_state: 随机种子
    
    返回:
        adata: AnnData 对象
    """
    np.random.seed(random_state)
    
    # 定义细胞类型（包含分支）
    cell_types = [
        "SSC",
        "Undiff.SPG",
        "Diff.SPG",
        "preL",
        "L",
        "Z",             # 分支点
        "L2",
        "L3",
        "preP",
        "P",
        "D",
        "SPC"
    ]
    
    n_cell_types = len(cell_types)
    cells_per_type = n_cells // n_cell_types
    
    # 生成 Y 形轨迹
    all_cells = []
    all_labels = []
    
    for i, ct in enumerate(cell_types):
        if ct in ["L", "Z"]:
            # 分支点：两个细胞类型在同一位置
            t = 0.5
            base_x = 0
            base_y = 0
        elif ct in ["L2", "L3", "preP", "P", "D", "SPC"]:
            # L 分支：向右上
            idx = cell_types.index(ct) - cell_types.index("L2")
            t = idx / 6 * np.pi / 3
            base_x = np.cos(t) * 4 + 2
            base_y = np.sin(t) * 4 + 2
        else:
            # 主干：向左下
            idx = cell_types.index(ct)
            t = idx / 4 * np.pi / 2
            base_x = -np.cos(t) * 5
            base_y = -np.sin(t) * 5
        
        cells_x = base_x + np.random.randn(cells_per_type) * 0.6
        cells_y = base_y + np.random.randn(cells_per_type) * 0.6
        
        all_cells.append(np.column_stack([cells_x, cells_y]))
        all_labels.extend([ct] * cells_per_type)
    
    X_umap = np.vstack(all_cells)
    
    adata = AnnData(
        X=np.random.randn(n_cells, 100),
        obs={'celltype': all_labels},
        obsm={'X_umap': X_umap}
    )
    
    # 定义带分支的树结构
    predefined_tree = {
        "SSC": ["Undiff.SPG"],
        "Undiff.SPG": ["Diff.SPG"],
        "Diff.SPG": ["preL"],
        "preL": ["L"],
        "L": ["Z", "L2"],  # 分支点
        "Z": ["preP"],     # Z 分支
        "L2": ["L3"],      # L 分支
        "L3": ["preP"],    # 汇合点
        "preP": ["P"],
        "P": ["D"],
        "D": ["SPC"]
    }
    
    adata.uns['predefined_tree'] = predefined_tree
    
    return adata


def plot_slingshot_result(slingshot_obj, adata, save_path=None, show=True):
    """
    绘制 Slingshot 结果
    
    参数:
        slingshot_obj: Slingshot 对象
        adata: 输入的 AnnData 数据
        save_path: 保存路径（可选）
        show: 是否显示（可选）
    
    返回:
        fig, axes: matplotlib 图形和轴
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    import seaborn as sns
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    
    # 1. 细胞类型分布
    celltype_key = 'celltype'
    unique_types = adata.obs[celltype_key].unique()
    n_types = len(unique_types)
    colors = sns.color_palette('tab20', n_types)
    color_map = {ct: colors[i] for i, ct in enumerate(unique_types)}
    
    for i, ct in enumerate(unique_types):
        mask = adata.obs[celltype_key] == ct
        axes[0, 0].scatter(
            adata.obsm['X_umap'][mask, 0],
            adata.obsm['X_umap'][mask, 1],
            c=[color_map[ct]],
            label=ct,
            s=20,
            alpha=0.7
        )
    
    axes[0, 0].set_title('Cell Types (Input)', fontsize=14, fontweight='bold')
    axes[0, 0].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    axes[0, 0].set_xlabel('UMAP 1')
    axes[0, 0].set_ylabel('UMAP 2')
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. 聚类中心
    sling = slingshot_obj
    cluster_colors = sns.color_palette('tab20', sling.num_clusters)
    
    axes[0, 1].scatter(
        sling.data[:, 0],
        sling.data[:, 1],
        c='lightgray',
        s=10,
        alpha=0.5,
        label='Cells'
    )
    axes[0, 1].scatter(
        sling.cluster_centres[:, 0],
        sling.cluster_centres[:, 1],
        c=cluster_colors,
        s=200,
        marker='X',
        edgecolors='black',
        linewidths=2,
        label='Cluster Centers'
    )
    
    # 标注起始节点
    axes[0, 1].scatter(
        sling.cluster_centres[sling.start_node, 0],
        sling.cluster_centres[sling.start_node, 1],
        c='red',
        s=300,
        marker='*',
        edgecolors='black',
        linewidths=2,
        label=f'Start Node ({sling.start_node})'
    )
    
    axes[0, 1].set_title('Cluster Centers & Start Node', fontsize=14, fontweight='bold')
    axes[0, 1].legend()
    axes[0, 1].set_xlabel('UMAP 1')
    axes[0, 1].set_ylabel('UMAP 2')
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. 树结构（MST 或预定义树）
    axes[0, 1].scatter(
        sling.data[:, 0],
        sling.data[:, 1],
        c='lightgray',
        s=10,
        alpha=0.3
    )
    
    # 绘制树边
    if hasattr(sling, 'tree') and sling.tree is not None:
        for parent, children in sling.tree.items():
            for child in children:
                axes[0, 1].plot(
                    [sling.cluster_centres[parent, 0], sling.cluster_centres[child, 0]],
                    [sling.cluster_centres[parent, 1], sling.cluster_centres[child, 1]],
                    'k-',
                    linewidth=2,
                    alpha=0.7
                )
    
    axes[0, 1].set_title('Tree Structure (MST or Predefined)', fontsize=14, fontweight='bold')
    
    # 4. 主曲线和谱系
    if sling.curves is not None:
        axes[1, 0].scatter(
            sling.data[:, 0],
            sling.data[:, 1],
            c='lightgray',
            s=10,
            alpha=0.5
        )
        
        # 绘制主曲线
        for l_idx, curve in enumerate(sling.curves):
            s_interp, p_interp, order = curve.unpack_params()
            axes[1, 0].plot(
                p_interp[order, 0],
                p_interp[order, 1],
                linewidth=3,
                label=f'Lineage {l_idx}',
                alpha=0.8
            )
        
        axes[1, 0].set_title('Principal Curves (Lineages)', fontsize=14, fontweight='bold')
        axes[1, 0].legend()
        axes[1, 0].set_xlabel('UMAP 1')
        axes[1, 0].set_ylabel('UMAP 2')
        axes[1, 0].grid(True, alpha=0.3)
    else:
        axes[1, 0].text(0.5, 0.5, 'Curves not fitted yet\nCall slingshot.fit()',
                       ha='center', va='center', fontsize=14)
        axes[1, 0].set_title('Principal Curves (Not Fitted)', fontsize=14, fontweight='bold')
    
    # 5. 伪时间分布
    if sling.curves is not None and hasattr(sling.curves[0], 'pseudotimes_interp'):
        # 计算统一伪时间
        pseudotime = np.zeros(sling.data.shape[0])
        
        for l_idx, lineage in enumerate(sling.lineages):
            curve = sling.curves[l_idx]
            cell_mask = np.logical_or.reduce(
                np.array([sling.cluster_label_indices == k for k in lineage])
            )
            pseudotime[cell_mask] = curve.pseudotimes_interp[cell_mask]
        
        # 归一化到 [0, 1]
        pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())
        
        scatter = axes[1, 1].scatter(
            sling.data[:, 0],
            sling.data[:, 1],
            c=pseudotime,
            cmap='viridis',
            s=20,
            alpha=0.7
        )
        
        plt.colorbar(scatter, ax=axes[1, 1], label='Pseudotime (normalized)')
        axes[1, 1].set_title('Pseudotime Distribution', fontsize=14, fontweight='bold')
        axes[1, 1].set_xlabel('UMAP 1')
        axes[1, 1].set_ylabel('UMAP 2')
        axes[1, 1].grid(True, alpha=0.3)
    else:
        axes[1, 1].text(0.5, 0.5, 'Pseudotime not computed yet\nCall slingshot.fit()',
                       ha='center', va='center', fontsize=14)
        axes[1, 1].set_title('Pseudotime (Not Computed)', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    if show:
        plt.show()
    
    return fig, axes


def plot_lineage_network(slingshot_obj, save_path=None, show=True):
    """
    绘制谱系网络图（使用 networkx）
    
    参数:
        slingshot_obj: Slingshot 对象
        save_path: 保存路径（可选）
        show: 是否显示（可选）
    
    返回:
        fig, ax: matplotlib 图形和轴
    """
    try:
        import networkx as nx
    except ImportError:
        print("networkx not installed. Install with: pip install networkx")
        return None, None
    
    fig, ax = plt.subplots(figsize=(12, 16))
    
    G = nx.DiGraph()
    
    # 添加节点和边
    if hasattr(slingshot_obj, 'tree') and slingshot_obj.tree is not None:
        for parent, children in slingshot_obj.tree.items():
            G.add_node(parent)
            for child in children:
                G.add_edge(parent, child)
    
    # 使用 hierarchical 布局
    pos = nx.drawing.nx_agraph.graphviz_layout(G, prog='dot')
    
    # 绘制节点
    node_sizes = [1000] * len(G.nodes())
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color='lightblue',
                          alpha=0.8, ax=ax)
    
    # 绘制边
    nx.draw_networkx_edges(G, pos, arrows=True, arrowsize=20,
                          edge_color='gray', width=2, ax=ax)
    
    # 绘制标签
    labels = {}
    for node in G.nodes():
        if hasattr(slingshot_obj, 'string_labels') and slingshot_obj.string_labels:
            labels[node] = slingshot_obj.int_to_celltype.get(node, str(node))
        else:
            labels[node] = str(node)
    
    nx.draw_networkx_labels(G, pos, labels, font_size=10,
                           bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='lightgray'),
                           ax=ax)
    
    ax.set_title('Lineage Network', fontsize=16, fontweight='bold')
    ax.axis('off')
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    if show:
        plt.show()
    
    return fig, ax


def run_test_pipeline(test_name="spermatogenesis"):
    """
    运行完整的测试流程
    
    参数:
        test_name: 测试名称 ("spermatogenesis" 或 "bifurcation")
    """
    import matplotlib.pyplot as plt
    
    print(f"\n{'='*60}")
    print(f"Running Slingshot Test: {test_name}")
    print(f"{'='*60}\n")
    
    # 1. 生成测试数据
    if test_name == "spermatogenesis":
        adata = generate_test_data_spermatogenesis(n_cells=520)
    elif test_name == "bifurcation":
        adata = generate_test_data_bifurcation(n_cells=600)
    else:
        raise ValueError(f"Unknown test_name: {test_name}")
    
    print(f"Generated {adata.n_obs} cells with {adata.n_vars} genes")
    print(f"Cell types: {adata.obs['celltype'].unique().tolist()}")
    print(f"Predefined tree: {list(adata.uns['predefined_tree'].keys())}\n")
    
    # 2. 初始化 Slingshot
    try:
        from hic_basic.slingshot import Slingshot
    except ImportError:
        print("ERROR: hic_basic.slingshot not available")
        print("Make sure the module is installed: pip install -e /path/to/hic_basic")
        return
    
    predefined_tree = adata.uns['predefined_tree']
    
    slingshot_obj = Slingshot(
        data=adata,
        celltype_key='celltype',
        obsm_key='X_umap',
        predefined_tree=predefined_tree,
        start_node="SSC",
        debug_level='verbose'
    )
    
    print(f"Slingshot initialized")
    print(f"Number of clusters: {slingshot_obj.num_clusters}")
    print(f"String labels: {slingshot_obj.string_labels}")
    print(f"Start node: {slingshot_obj.start_node}")
    if slingshot_obj.string_labels:
        print(f"Cell type mapping: {slingshot_obj.int_to_celltype}\n")
    
    # 3. 查看树结构
    print("Tree structure:")
    if hasattr(slingshot_obj, 'label_tree') and slingshot_obj.label_tree is not None:
        for parent, children in slingshot_obj.label_tree.items():
            children_str = [slingshot_obj.int_to_celltype.get(c, str(c)) 
                          if slingshot_obj.string_labels else str(c) 
                          for c in children]
            print(f"  {parent} -> {children_str}")
    print()
    
    # 4. 拟合主曲线
    print("Fitting principal curves...")
    slingshot_obj.fit(num_epochs=10)
    print(f"Number of lineages: {len(slingshot_obj.lineages)}")
    print(f"Lineages: {[list(l) for l in slingshot_obj.lineages]}\n")
    
    # 5. 列出谱系
    print("Lineage details:")
    slingshot_obj.list_lineages()
    print()
    
    # 6. 绘图
    print("Generating plots...")
    
    # 主图
    plot_slingshot_result(
        slingshot_obj,
        adata,
        save_path=f'/tmp/slingshot_test_{test_name}_main.png',
        show=False
    )
    
    # 网络图
    plot_lineage_network(
        slingshot_obj,
        save_path=f'/tmp/slingshot_test_{test_name}_network.png',
        show=False
    )
    
    print(f"\nTest completed successfully!")
    print(f"Plots saved to /tmp/slingshot_test_{test_name}_*.png")
    print(f"{'='*60}\n")
    
    return slingshot_obj, adata


if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')  # 非交互式后端
    import matplotlib.pyplot as plt
    
    # 运行测试
    sling, adata = run_test_pipeline("spermatogenesis")
    # sling2, adata2 = run_test_pipeline("bifurcation")
