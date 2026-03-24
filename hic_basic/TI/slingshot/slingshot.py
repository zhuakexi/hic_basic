"""
Slingshot - 细胞谱系和伪时间推断

基于 pyslingshot 的增强版本，支持：
- 预定义树结构 (predefined_tree)
- 字符串标签 (string_labels)
- 缺失节点自动处理
"""

from typing import Union, Dict, List, Optional

import numpy as np
from anndata import AnnData

from pcurvepy2 import PrincipalCurve
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.interpolate import interp1d
from sklearn.neighbors import KernelDensity
from collections import deque
from tqdm.autonotebook import tqdm

from .utils import (
    mahalanobis,
    isint,
    isstr,
    scale_to_range,
    infer_cluster_label_indices
)
from .lineage import Lineage


class Slingshot:
    """
    Slingshot 细胞谱系推断类
    
    基于主曲线 (Principal Curve) 的轨迹推断算法，支持：
    - MST 自动构建谱系树
    - 预定义树结构（整合生物学先验知识）
    - 字符串细胞类型标签
    - 缺失细胞类型自动处理
    """
    
    def __init__(
            self,
            data: Union[AnnData, np.ndarray],
            cluster_labels_onehot=None,
            celltype_key: Optional[str] = None,
            obsm_key='X_umap',
            start_node: Union[str, int] = 0,
            end_nodes: Optional[List[Union[str, int]]] = None,
            predefined_tree: Optional[Dict[str, List[str]]] = None,
            debug_level=None
    ):
        """
        初始化 Slingshot 对象
        
        参数:
            data: AnnData 对象或 numpy 数组，包含降维后的数据 (num_cells, 2)
            cluster_labels_onehot: one-hot 编码的聚类标签，仅在 data 不是 AnnData 时需要
            celltype_key: AnnData.obs 中的细胞类型列名
            obsm_key: AnnData.obsm 中的降维数据键
            start_node: 起始节点（字符串或整数）
            end_nodes: 终端节点列表
            predefined_tree: 预定义的树结构（字典映射父细胞类型到子细胞类型列表）
            debug_level: 调试级别
        """
        self.predefined_tree = predefined_tree
        self.string_labels = False
        self.celltype_to_int = {}
        self.int_to_celltype = {}

        if isinstance(data, AnnData):
            assert celltype_key is not None, "Must provide celltype key if data is an AnnData object"
            cluster_labels = data.obs[celltype_key]

            # 检查标签是否为字符串并创建映射
            if cluster_labels.dtype == object or isinstance(cluster_labels.iloc[0], str):
                self.string_labels = True
                unique_labels = sorted(cluster_labels.unique())
                self.celltype_to_int = {label: i for i, label in enumerate(unique_labels)}
                self.int_to_celltype = {i: label for label, i in self.celltype_to_int.items()}
                cluster_labels_numeric = np.array([self.celltype_to_int[label] for label in cluster_labels])
            else:
                cluster_labels_numeric = cluster_labels.values

            self.cluster_label_indices = infer_cluster_label_indices(cluster_labels_numeric)
            cluster_labels_onehot = np.zeros((cluster_labels_numeric.shape[0], self.cluster_label_indices.max() + 1))
            cluster_labels_onehot[np.arange(cluster_labels_numeric.shape[0]), self.cluster_label_indices] = 1

            data = data.obsm[obsm_key]
        else:
            assert cluster_labels_onehot is not None, "Must provide cluster labels if data is not an AnnData object"
            cluster_labels = cluster_labels_onehot.argmax(axis=1)
            self.cluster_label_indices = infer_cluster_label_indices(cluster_labels)

        self.data = data
        self.cluster_labels_onehot = cluster_labels_onehot
        self.cluster_labels = cluster_labels
        self.num_clusters = self.cluster_label_indices.max() + 1

        # 获取数据中可用的细胞类型
        self.available_celltypes = set(range(self.num_clusters))
        if self.string_labels:
            self.available_celltype_strings = set(self.int_to_celltype.values())
        else:
            self.available_celltype_strings = None

        # 转换 start_node 和 end_nodes 为内部整数表示
        if self.string_labels:
            if isinstance(start_node, str):
                if start_node not in self.celltype_to_int:
                    raise ValueError(f"Start node '{start_node}' not found in data cell types")
                self.start_node = self.celltype_to_int[start_node]
            else:
                self.start_node = start_node

            if end_nodes is None:
                self.end_nodes = []
            else:
                self.end_nodes = []
                for node in end_nodes:
                    if isinstance(node, str):
                        if node in self.celltype_to_int:
                            self.end_nodes.append(self.celltype_to_int[node])
                        else:
                            print(f"Warning: End node '{node}' not found in data, skipping")
                    else:
                        self.end_nodes.append(node)
        else:
            self.start_node = start_node
            self.end_nodes = [] if end_nodes is None else end_nodes

        cluster_centres = [data[self.cluster_label_indices == k].mean(axis=0) for k in range(self.num_clusters)]
        self.cluster_centres = np.stack(cluster_centres)
        self.lineages = None      # list of Lineages
        self.cluster_lineages = None # lineages belonging to each cluster
        self.curves = None   # list of principle curves len = #lineages
        self.cell_weights = None  # weights indicating cluster assignments
        self.distances = None
        self.branch_clusters = None
        self._tree = None

        # 绘图和打印
        debug_level = 0 if debug_level is None else dict(verbose=1)[debug_level]
        self.debug_level = debug_level
        self._set_debug_axes(None)

        # 构建收缩步骤的平滑核
        self.kernel_x = np.linspace(-3, 3, 512)
        kde = KernelDensity(bandwidth=1., kernel='gaussian')
        kde.fit(np.zeros((self.kernel_x.shape[0], 1)))
        self.kernel_y = np.exp(kde.score_samples(self.kernel_x.reshape(-1, 1)))

    @property
    def tree(self):
        """获取树结构（MST 或预定义树）"""
        if self._tree is None:
            if self.predefined_tree is not None:
                self.construct_tree_from_predefined()
            else:
                self.construct_mst(self.start_node)
        return self._tree

    @property
    def label_tree(self):
        """返回带字符串标签的树结构"""
        if not hasattr(self, '_tree') or self._tree is None:
            return None

        if self.string_labels:
            string_tree = {}
            for int_parent, int_children in self._tree.items():
                parent_str = self.int_to_celltype[int_parent]
                children_str = [self.int_to_celltype[child] for child in int_children]
                string_tree[parent_str] = children_str
            return string_tree
        return self._tree

    def construct_tree_from_predefined(self):
        """从预定义树结构构建，处理缺失节点"""
        if self.string_labels:
            complete_tree = {}
            for parent, children in self.predefined_tree.items():
                complete_tree[parent] = children

            pruned_tree = {}
            visited = set()

            def build_pruned_tree(current_node):
                if current_node in visited:
                    return
                visited.add(current_node)

                if current_node not in self.celltype_to_int:
                    print(f"Warning: Cell type '{current_node}' in predefined tree not found in data, skipping")
                    if current_node in complete_tree:
                        for child in complete_tree[current_node]:
                            build_pruned_tree(child)
                    return

                current_int = self.celltype_to_int[current_node]
                if current_int not in pruned_tree:
                    pruned_tree[current_int] = []

                if current_node in complete_tree:
                    for child in complete_tree[current_node]:
                        if child in self.celltype_to_int:
                            child_int = self.celltype_to_int[child]
                            pruned_tree[current_int].append(child_int)
                            build_pruned_tree(child)
                        else:
                            print(f"Warning: Cell type '{child}' in predefined tree not found in data, skipping")
                            if child in complete_tree:
                                for grandchild in complete_tree[child]:
                                    if grandchild in self.celltype_to_int:
                                        grandchild_int = self.celltype_to_int[grandchild]
                                        pruned_tree[current_int].append(grandchild_int)
                                        build_pruned_tree(grandchild)

            start_node_str = self.int_to_celltype[self.start_node]
            build_pruned_tree(start_node_str)
            self._tree = pruned_tree

        else:
            complete_tree = {}
            for parent, children in self.predefined_tree.items():
                parent_int = int(parent) if isinstance(parent, str) else parent
                children_int = [int(child) if isinstance(child, str) else child for child in children]
                complete_tree[parent_int] = children_int

            pruned_tree = {}
            visited = set()

            def build_pruned_tree(current_node):
                if current_node in visited:
                    return
                visited.add(current_node)

                if current_node not in self.available_celltypes:
                    print(f"Warning: Cell type {current_node} in predefined tree not found in data, skipping")
                    if current_node in complete_tree:
                        for child in complete_tree[current_node]:
                            build_pruned_tree(child)
                    return

                if current_node not in pruned_tree:
                    pruned_tree[current_node] = []

                if current_node in complete_tree:
                    for child in complete_tree[current_node]:
                        if child in self.available_celltypes:
                            pruned_tree[current_node].append(child)
                            build_pruned_tree(child)
                        else:
                            print(f"Warning: Cell type {child} in predefined tree not found in data, skipping")
                            if child in complete_tree:
                                for grandchild in complete_tree[child]:
                                    if grandchild in self.available_celltypes:
                                        pruned_tree[current_node].append(grandchild)
                                        build_pruned_tree(grandchild)

            build_pruned_tree(self.start_node)
            self._tree = pruned_tree

        if self.start_node not in self._tree:
            raise ValueError(f"Start node {self.start_node} not found in the constructed tree")

        if self.debug_level > 0:
            print("Constructed tree structure:")
            for parent, children in self._tree.items():
                parent_str = self.int_to_celltype[parent] if self.string_labels else str(parent)
                children_str = [self.int_to_celltype[child] if self.string_labels else str(child) for child in children]
                print(f"  {parent_str} -> {children_str}")

    def load_params(self, filepath):
        """加载参数"""
        if self.curves is None:
            self.get_lineages()
        params = np.load(filepath, allow_pickle=True).item()
        self.curves = params['curves']
        self.cell_weights = params['cell_weights']
        self.distances = params['distances']

    def save_params(self, filepath):
        """保存参数"""
        params = dict(
            curves=self.curves,
            cell_weights=self.cell_weights,
            distances=self.distances
        )
        np.save(filepath, params)

    def _set_debug_axes(self, axes):
        self.debug_axes = axes
        self.debug_plot_mst = axes is not None
        self.debug_plot_lineages = axes is not None
        self.debug_plot_avg = axes is not None

    def construct_mst(self, start_node):
        """
        构建最小生成树
        
        参数:
            start_node: 起始节点
        返回:
            children: 字典映射聚类到子节点
        """
        emp_covs = np.stack([np.cov(self.data[self.cluster_label_indices == i].T) for i in range(self.num_clusters)])
        dists = np.zeros((self.num_clusters, self.num_clusters))
        for i in range(self.num_clusters):
            for j in range(i, self.num_clusters):
                dist = mahalanobis(
                    self.cluster_centres[i],
                    self.cluster_centres[j],
                    emp_covs[i],
                    emp_covs[j]
                )
                dists[i, j] = dist
                dists[j, i] = dist

        mst_dists = np.delete(np.delete(dists, self.end_nodes, axis=0), self.end_nodes, axis=1)
        tree = minimum_spanning_tree(mst_dists)
        index_mapping = np.array([c for c in range(self.num_clusters - len(self.end_nodes))])
        for i, end_node in enumerate(self.end_nodes):
            index_mapping[end_node - i:] += 1

        connections = {k: list() for k in range(self.num_clusters)}
        cx = tree.tocoo()
        for i, j, v in zip(cx.row, cx.col, cx.data):
            i = index_mapping[i]
            j = index_mapping[j]
            connections[i].append(j)
            connections[j].append(i)

        for end in self.end_nodes:
            i = np.argmin(np.delete(dists[end], self.end_nodes))
            connections[i].append(end)
            connections[end].append(i)

        visited = [False for _ in range(self.num_clusters)]
        queue = list()
        queue.append(start_node)
        children = {k: list() for k in range(self.num_clusters)}
        while len(queue) > 0:
            current_node = queue.pop()
            visited[current_node] = True
            for child in connections[current_node]:
                if not visited[child]:
                    children[current_node].append(child)
                    queue.append(child)

        if self.debug_plot_mst:
            # Debug plotting disabled - use plot_slingshot() for visualization
            self.debug_plot_mst = False

        self._tree = children
        return children

    def get_lineages(self):
        """从 MST 提取谱系"""
        tree = self.construct_mst(self.start_node)

        branch_clusters = deque()
        def recurse_branches(path, v):
            num_children = len(tree[v])
            if num_children == 0:
                return path + [v, None]
            elif num_children == 1:
                return recurse_branches(path + [v], tree[v][0])
            else:
                branch_clusters.append(v)
                return [recurse_branches(path + [v], tree[v][i]) for i in range(num_children)]

        def flatten(li):
            if li[-1] is None:
                yield Lineage(li[:-1])
            else:
                for l in li:
                    yield from flatten(l)

        lineages = recurse_branches([], self.start_node)
        lineages = list(flatten(lineages))
        self.lineages = lineages
        self.branch_clusters = branch_clusters

        self.cluster_lineages = {k: list() for k in range(self.num_clusters)}
        for l_idx, lineage in enumerate(self.lineages):
            for k in lineage:
                self.cluster_lineages[k].append(l_idx)

        if self.debug_level > 0:
            print('Lineages:', lineages)

    def fit(self, num_epochs=10, debug_axes=None):
        """
        拟合主曲线
        
        参数:
            num_epochs: 迭代次数
            debug_axes: 调试绘图轴
        """
        self._set_debug_axes(debug_axes)
        if self.curves is None:
            self.get_lineages()
            self.construct_initial_curves()
            self.cell_weights = [self.cluster_labels_onehot[:, self.lineages[l].clusters].sum(axis=1)
                                 for l in range(len(self.lineages))]
            self.cell_weights = np.stack(self.cell_weights, axis=1)

        for epoch in tqdm(range(num_epochs)):
            self.calculate_cell_weights()
            self.fit_lineage_curves()

            for l_idx, lineage in enumerate(self.lineages):
                curve = self.curves[l_idx]
                min_time = np.min(curve.pseudotimes_interp[self.cell_weights[:, l_idx] > 0])
                curve.pseudotimes_interp -= min_time

            shrinkage_percentages, cluster_children, cluster_avg_curves = self.avg_curves()
            self.shrink_curves(cluster_children, shrinkage_percentages, cluster_avg_curves)

            self.debug_plot_lineages = False
            self.debug_plot_avg = False

            # Debug plotting disabled - use plot_slingshot() for visualization

    def construct_initial_curves(self):
        """使用分段线性初始化构建谱系主曲线"""
        piecewise_linear = list()
        distances = list()

        for l_idx, lineage in enumerate(self.lineages):
            p = np.stack(self.cluster_centres[lineage.clusters])
            s = np.zeros(p.shape[0])

            cell_mask = np.logical_or.reduce(
                np.array([self.cluster_label_indices == k for k in lineage]))
            cells_involved = self.data[cell_mask]

            curve = PrincipalCurve(k=3)
            curve.project_to_curve(cells_involved, points=p)
            d_sq, dist = curve.project_to_curve(self.data, points=curve.points_interp[curve.order])
            distances.append(d_sq)

            piecewise_linear.append(curve)

        self.curves = piecewise_linear
        self.distances = distances

    def fit_lineage_curves(self):
        """使用三次样条和数据投影更新曲线"""
        assert self.lineages is not None
        assert self.curves is not None
        distances = list()

        for l_idx, lineage in enumerate(self.lineages):
            curve = self.curves[l_idx]
            curve.fit(self.data, max_iter=1, w=self.cell_weights[:, l_idx])

            if self.debug_plot_lineages:
                cell_mask = np.logical_or.reduce(
                    np.array([self.cluster_label_indices == k for k in lineage]))
                cells_involved = self.data[cell_mask]
                self.debug_axes[0, 1].scatter(cells_involved[:, 0], cells_involved[:, 1], s=2, alpha=0.5)
                alphas = curve.pseudotimes_interp
                alphas = (alphas - alphas.min()) / (alphas.max() - alphas.min())
                for i in np.random.permutation(self.data.shape[0])[:50]:
                    path_from = (self.data[i][0], curve.points_interp[i][0])
                    path_to = (self.data[i][1], curve.points_interp[i][1])
                    self.debug_axes[0, 1].plot(path_from, path_to, c='black', alpha=alphas[i])
                self.debug_axes[0, 1].plot(curve.points_interp[curve.order, 0],
                                           curve.points_interp[curve.order, 1], label=str(lineage))

            d_sq, dist = curve.project_to_curve(self.data, curve.points_interp[curve.order])
            distances.append(d_sq)
        self.distances = distances
        if self.debug_plot_lineages:
            self.debug_axes[0, 1].legend()

    def calculate_cell_weights(self):
        """计算细胞权重"""
        cell_weights = [self.cluster_labels_onehot[:, self.lineages[l].clusters].sum(axis=1)
                        for l in range(len(self.lineages))]
        cell_weights = np.stack(cell_weights, axis=1)

        d_sq = np.stack(self.distances, axis=1)
        d_ord = np.argsort(d_sq, axis=None)
        w_prob = cell_weights/cell_weights.sum(axis=1, keepdims=True)
        w_rnk_d = np.cumsum(w_prob.reshape(-1)[d_ord]) / w_prob.sum()

        z = d_sq
        z_shape = z.shape
        z = z.reshape(-1)
        z[d_ord] = w_rnk_d
        z = z.reshape(z_shape)
        z_prime = 1 - z ** 2
        z_prime[cell_weights == 0] = np.nan
        w0 = cell_weights.copy()
        cell_weights = z_prime / np.nanmax(z_prime, axis=1, keepdims=True)
        np.nan_to_num(cell_weights, nan=1, copy=False)
        cell_weights[cell_weights > 1] = 1
        cell_weights[cell_weights < 0] = 0
        cell_weights[w0 == 0] = 0

        reassign = True
        if reassign:
            cell_weights[z < .5] = 1
            ridx = (z.max(axis=1) > .9) & (cell_weights.min(axis=1) < .1)
            w0 = cell_weights[ridx]
            z0 = z[ridx]
            w0[(z0 > .9) & (w0 < .1)] = 0
            cell_weights[ridx] = w0

        self.cell_weights = cell_weights

    def avg_curves(self):
        """从叶子到根计算分支的平均曲线"""
        cell_weights = self.cell_weights
        shrinkage_percentages = list()
        cluster_children = dict()
        lineage_avg_curves = dict()
        cluster_avg_curves = dict()
        branch_clusters = self.branch_clusters.copy()
        
        if self.debug_level > 0:
            print('Reversing from leaf to root')
        # Debug plotting disabled - use plot_slingshot() for visualization

        while len(branch_clusters) > 0:
            k = branch_clusters.pop()
            branch_lineages = self.cluster_lineages[k]
            cluster_children[k] = set()
            for l_idx in branch_lineages:
                if l_idx in lineage_avg_curves:
                    curve = lineage_avg_curves[l_idx]
                else:
                    curve = self.curves[l_idx]
                cluster_children[k].add(curve)

            branch_curves = list(cluster_children[k])
            if self.debug_level > 0:
                print(f'Averaging branch @{k} with lineages:', branch_lineages, branch_curves)

            avg_curve = self.avg_branch_curves(branch_curves)
            cluster_avg_curves[k] = avg_curve

            common = cell_weights[:, branch_lineages] > 0
            common_mask = common.mean(axis=1) == 1.
            shrinkage_percent = dict()
            for curve in branch_curves:
                shrinkage_percent[curve] = self.shrinkage_percent(curve, common_mask)
            shrinkage_percentages.append(shrinkage_percent)

            for l in branch_lineages:
                lineage_avg_curves[l] = avg_curve

        if self.debug_plot_avg:
            self.debug_axes[1, 0].legend()
        return shrinkage_percentages, cluster_children, cluster_avg_curves

    def shrink_curves(self, cluster_children, shrinkage_percentages, cluster_avg_curves):
        """从根到叶收缩曲线"""
        branch_clusters = self.branch_clusters.copy()
        while len(branch_clusters) > 0:
            k = branch_clusters.popleft()
            shrinkage_percent = shrinkage_percentages.pop()
            branch_curves = list(cluster_children[k])
            cluster_avg_curve = cluster_avg_curves[k]
            if self.debug_level > 0:
                print(f'Shrinking branch @{k} with curves:', branch_curves)
            self.shrink_branch_curves(branch_curves, cluster_avg_curve, shrinkage_percent)

    def shrink_branch_curves(self, branch_curves, avg_curve, shrinkage_percent):
        """收缩分支曲线到平均曲线"""
        num_dims_reduced = branch_curves[0].points_interp.shape[1]

        for curve in branch_curves:
            pct = shrinkage_percent[curve]
            s_interp, p_interp, order = curve.unpack_params()
            avg_s_interp, avg_p_interp, avg_order = avg_curve.unpack_params()
            shrunk_curve = np.zeros_like(p_interp)
            for j in range(num_dims_reduced):
                orig = p_interp[order, j]
                avg = np.interp(s_interp[order], avg_s_interp[avg_order], avg_p_interp[avg_order, j])
                shrunk_curve[:, j] = (avg * pct + orig * (1 - pct))
            curve.project_to_curve(self.data, points=shrunk_curve)

    def shrinkage_percent(self, curve, common_ind):
        """确定曲线收缩百分比"""
        s_interp, order = curve.pseudotimes_interp, curve.order
        x = self.kernel_x
        y = self.kernel_y
        y = (y.sum() - np.cumsum(y)) / sum(y)
        q1 = np.percentile(s_interp[common_ind], 25)
        q3 = np.percentile(s_interp[common_ind], 75)
        a = q1 - 1.5 * (q3 - q1)
        b = q3 + 1.5 * (q3 - q1)
        x = scale_to_range(x, a=a, b=b)
        if q1 == q3:
            pct_l = np.zeros(s_interp.shape[0])
        else:
            pct_l = np.interp(s_interp[order], x, y)
        return pct_l

    def avg_branch_curves(self, branch_curves):
        """计算分支曲线的平均曲线"""
        num_cells = branch_curves[0].points_interp.shape[0]
        num_dims_reduced = branch_curves[0].points_interp.shape[1]

        branch_s_interps = np.stack([c.pseudotimes_interp for c in branch_curves], axis=1)
        max_shared_pseudotime = branch_s_interps.max(axis=0).min()
        combined_pseudotime = np.linspace(0, max_shared_pseudotime, num_cells)
        curves_dense = list()
        
        for curve in branch_curves:
            lineage_curve = np.zeros((combined_pseudotime.shape[0], num_dims_reduced))
            order = curve.order
            for j in range(num_dims_reduced):
                lin_interpolator = interp1d(
                    curve.pseudotimes_interp[order],
                    curve.points_interp[order, j],
                    assume_sorted=True
                )
                lineage_curve[:, j] = lin_interpolator(combined_pseudotime)
            curves_dense.append(lineage_curve)

        curves_dense = np.stack(curves_dense, axis=1)
        avg = curves_dense.mean(axis=1)
        
        avg_curve = PrincipalCurve()
        avg_curve.project_to_curve(self.data, points=avg)
        
        if self.debug_plot_avg:
            self.debug_axes[1, 0].plot(avg[:, 0], avg[:, 1], c='blue', linestyle='--', label='average', alpha=0.7)
            _, p_interp, order = avg_curve.unpack_params()
            self.debug_axes[1, 0].plot(p_interp[order, 0], p_interp[order, 1], c='red', label='data projected', alpha=0.7)

        return avg_curve

    def list_lineages(self, cluster_to_label=None):
        """列出谱系（支持字符串标签）"""
        if cluster_to_label is None and self.string_labels:
            cluster_to_label = self.int_to_celltype

        for lineage in self.lineages:
            if cluster_to_label is not None:
                lineage_labels = [cluster_to_label[l] for l in lineage]
                print(', '.join(lineage_labels))
            else:
                print(', '.join([str(l) for l in lineage]))

    @property
    def unified_pseudotime(self):
        """获取统一的伪时间"""
        pseudotime = np.zeros_like(self.curves[0].pseudotimes_interp)
        for l_idx, lineage in enumerate(self.lineages):
            curve = self.curves[l_idx]
            cell_mask = np.logical_or.reduce(
                np.array([self.cluster_label_indices == k for k in lineage])
            )
            pseudotime[cell_mask] = curve.pseudotimes_interp[cell_mask]
        return pseudotime
