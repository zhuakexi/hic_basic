# --- meta data: group ---
from functools import total_ordering
from io import StringIO
import re
import yaml

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def drange(a,b,seats):
    out = []
    for i in range(a,b):
        fill_0 = seats - len(str(i))
        if fill_0 < 0:
            out.append(str(i))
        else:
            out.append("0"*fill_0 + str(i))
    return out
def parse_group_string(group_string, last_mouse):
    """
    Parsing the cell meta data record from experiment.
    Input:
        group_string: string to be parsed, should be in yaml format.
        last_mouse: largest mouse number in past experiments 
    """
    groups = []
    
    f = StringIO(group_string)
    kv = yaml.load(f, yaml.BaseLoader)
    f.close()
    
    for group in kv.keys():
        # group: samples with same prefix
        current_cell_prefix = group[0:8]
        current_digits = len(group.split("-")[1])
        for series in kv[group].keys():
            # series: samples with same group_name
            start = int(series.split("-")[0])
            end = int(series.split("-")[1]) + 1
            partition = kv[group][series]["p"]
            cell_type = kv[group][series]["c"]
            mouse_number = int(kv[group][series]["m"])
            clock = kv[group][series]["o"]
            
            for i in drange(start, end, current_digits):
                # sample_name, group, partition, cell_type
                i_sample_name = current_cell_prefix + i
                i_group = "m" + str(last_mouse + mouse_number) + "_c" + cell_type + "_o" + clock
                i_partition = partition
                i_cell_type = "c" + cell_type
                groups.append((i_sample_name, i_group, i_partition, i_cell_type))
    groups = pd.DataFrame(groups, columns = ["sample_name","group","partition","cell_type"])
    groups = groups.set_index("sample_name")
    if not groups.index.is_unique:
        raise ValueError("parse_group_string Error: dranges overlap.")
    return groups

# --- meta data: group: parser ---
@total_ordering
class ExpTime:
    norm_grp_mapper={
    'm9_c1_o23':"m9_c1_o2300_d0",
    'm12_c1_o2130':"m12_c1_o2130_d0",
    'm46_c1orc2_o20302100':"m46_c1orc2_o20302100_d0",
    'm46_c1_o20302100':"m46_c1_o20302100_d0",
    'm50_c1_o16001700':"m50_c1_o16001700_d0",
    'm51_c1_o16001700':"m51_c1_o16001700_d0",
    'm53_c1_o21452200':"m53_c1_o21452200_d0",
    "m22_c2_o20302100": "m22_c2_o20302100_d0",
    "m12_c2_o2130": "m12_c2_o2130_d0",
    "m13_c2_o23": "m13_c2_o2300_d0",
    "m23_c2_o22002230": "m23_c2_o22002230_d0",
    "m9_c2_o23": "m9_c2_o2300_d0",
    "m27_c2_o700730": "m27_c2_o07000730_d1",
    "m10_c2_o1": "m10_c2_o0100_d1",
    "m20_c2_o530": "m20_c2_o0530_d1",
    "m4_c2_o1112": "m4_c2_o11001200_d1",
    "m18_c2_o500530": "m18_c2_o05000530_d1",
    "m29_c2_o730": "m29_c2_o0730_d1",
    "m21_c2_o530": "m21_c2_o0530_d1",
    "m19_c2_o630700": "m19_c2_o06300700_d1",
    "m1r_c2_o11": "m1r_c2_o1100_d1",
    "m24_c2_o10": "m24_c2_o1000_d1",
    "m28_c2_o930": "m28_c2_o0930_d1",
    "m30_c2_o1330": "m30_c2_o1330_d1",
    "m31_c2_o1011": "m31_c2_o10001100_d1",
    "m26_c2_o12": "m26_c2_o1200_d1",
    "m7_c2_o10": "m7_c2_o1000_d1",
    "m1l_c2_o10": "m1l_c2_o1000_d1",
    "m3_c2_o1011": "m3_c2_o10001100_d1",
    "m2_c2_o16": "m2_c2_o1600_d1",
    "m32_c2_o15301600" : "m32_c2_o15301600_d1",
    "m33_c2_o15301600" : "m33_c2_o15301600_d1",
    "m34_c2_o1430" : "m34_c2_o1430_d1",
    "m35_c2_o15301600" : "m35_c2_o15301600_d1",
    "c2_m36_o21002130" : "c2_m36_o21002130_d0",
    "c2_m37_o02000300" : "c2_m37_o02000300_d1",
    "m46_c2_o20302100" : "m46_c2_o20302100_d0",
    "m46_c1orc2_o20302100" : "m46_c1orc2_o20302100_d0",
    "m52_c2_o20152045" : "m52_c2_o20152045_d0",
    "m53_c2_o21452200" : "m53_c2_o21452200_d0",
    "m47_c2_o20002030" : "m47_c2_o20002030_d0",
    "m48_c2_o21002130" : "m48_c2_o21002130_d0",
    "m6_c4_o14" : "m6_c4_o1400_d1",
    "m2_c4_o16" : "m2_c4_o1600_d1",
    "m2_c4_o19" : "m2_c4_o1900_d1",
    "m5_c4_o1112" : "m5_c4_o11001200_d1",
    "m15_c4_o20" : "m15_c4_o2000_d1",
    "m16_c4_o2223" : "m16_c4_o22002300_d1",
    "m25_c4_o21p3" : "m25_c4_o2100p3_d1",
    "m26_c4_o12" : "m26_c4_o1200_d1",
    "m35_c4_o15301600" : "m35_c4_o15301600_d1",
    "c4_m38_d3o04000400" : "c4_m38_o04000400_d2",
    "m41_c4_o22002230" : "m41_c4_o22002230_d1",
    "m42_c4_o00000030" : "m42_c4_o00000030_d2",
    "m40_c4_o0800830" : "m40_c4_o08000830_d2",
    'm33_c3or4_o15301600':'m33_c3or4_o15301600_d1',
    'c8to16_m39_d3o08300830':'c8to16_m39_o08300830_d2',
    'c8to16_m38_d3o04000400':'c8to16_m38_o04000400_d2',
    'm40_c8_o08000830':'m40_c8_o08000830_d2',
    "m49_c8_o10301030":"m49_c8_o10301030_d2"
    }
    @staticmethod
    def d2(time_str):
        if len(time_str) != 2:
            raise ValueError("d2: length must be 2")
        if time_str == "00":
            return 0
        return int(time_str.lstrip("0"))
    @staticmethod
    def find_time_str(group):
        """
        Find time substring in sample's group_string.
        """
        ts_re = re.compile("[_,^]+o(\d+)([a-z,A-Z]\d)*[$,_]+")
        # search only returns first matching.
        ts = ts_re.search(group)
        if ts is None:
            raise ValueError("Can't find time string in group name" + group)
        else:
            if ts.group(2) is not None:
                print("[find_time_str] Note: {} has additional clock annotation: {}.".format(group,ts.group(2)))
            return ts.group(1)
    @staticmethod
    def find_day(group):
        """
        Find day mark.
        """
        day_re = re.compile("d(\d+)")
        day = day_re.search(group)
        if day is None:
            raise ValueError("Can't find day mark in group name " + group)
        else:
            return int(day.group(1))
    def __init__(self, group):
        if group in self.norm_grp_mapper:
            norm_grp = self.norm_grp_mapper[group]
        else:
            norm_grp = group
        time_str = self.find_time_str(norm_grp)
        self.day = self.find_day(norm_grp)
        if len(time_str) == 8:
            self.tah = self.d2(time_str[0:2])
            self.tam = self.d2(time_str[2:4])
            self.tbh = self.d2(time_str[4:6])
            self.tbm = self.d2(time_str[6:8])
        elif len(time_str) == 4:
            self.tah = self.d2(time_str[0:2])
            self.tam = self.d2(time_str[2:4])
            self.tbh = self.tah
            self.tbm = self.tam
        else:
            print(time_str)
            raise ValueError("time_str must be either 8 or 4 chrs: " + group)
        self.hpf = (self.day*24*60 + self.tah*60 + self.tam)/60
    def __eq__(self, other):
        return (self.tah == other.tah) and (self.tam == other.tam) \
            and (self.tbh == other.tbh) and (self.tbm == other.tbm) \
            and (self.day == other.day)
    def __lt__(self, other):
        #return ((self.tah < other.tah) and (self.tam < other.tam)) \
        #    or ((self.tbh < other.tbh) and (self.tbm < other.tbm))
        #return self.tah < other.tah and self.tam < other.tam
        return (self.day*24*60 + self.tah*60 + self.tam) < (other.day*24*60 + other.tah*60 + other.tam)
    def __str__(self):
        return "day{} {}:{} {}:{}".format(self.day, self.tah, self.tam, self.tbh, self.tbm)
def add_group_order(annote):
    """
    Add additional order index col indicating sample collection time order to input df.
    Input must have a "group" col.
    """
    sorted_grp = sorted(annote["group"].unique(),key=ExpTime)
    group_orderi = pd.Series(
        list(range(len(sorted_grp))),
        index=sorted_grp,
        name="group_order_index"
    )
    if "group_order_index" in annote.columns:
        annote = annote.drop("group_order_index",axis=1)
    return pd.merge(annote,group_orderi,left_on="group",right_index=True)
def add_group_hour(annote, col="collect_hour"):
    """
    Add additional order index col to input df indicating hour post fertilizing time when the sample was collected.
    Input:
        annote: sample-feature table, must have a "group" col.
        col: name of new column.
    """
    grps = annote["group"].unique()
    grp_hour = pd.Series([ExpTime(grp).hpf for grp in grps], index=grps, name=col)
    if col in annote.columns:
        annote = annote.drop(col, axis=1)
    return pd.merge(annote, grp_hour, left_on="group",right_index=True)
def add_cell_type(annote):
    """
    Add cell_type col to input df according to group name.
    Input must have "group" col.
    """
    def get_cell_type(grp):
        cell_type_re = re.compile('(c[a-z,0-9]+)')
        cell_type = cell_type_re.search(grp)
        if cell_type is None:
            raise ValueError("[get_cell_type]: can't find cell_type mark in " + grp)
        else:
            return cell_type.group(1)
    if "cell_type" in annote.columns:
        annote = annote.drop("cell_type",axis=1)
    annote = annote.assign(cell_type = [get_cell_type(i) for i in annote["group"]])
    return annote

def plot_real_pseudo_time(adata, ps_col="velocity_pseudotime",hour=True):
    """
    Plot experiment time with pseudotime.
    (adata plot; scatter plot;)
    Input:
        adata: AnnData object with "group" col in obs.
        ps_col: obs col storing pseudotime to plot real time with. 
        hour: sort by real experiment hour; False to use relative ordering.
    """
    adata.obs = add_group_hour(adata.obs, "collect_hour")
    data = adata.obs.sort_values(ps_col)
    fig = px.scatter(
        x = data.index,
        y = data["collect_hour"],
        color = data["group"]
    )
    return fig