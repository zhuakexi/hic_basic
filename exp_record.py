# --- meta data: group ---
from functools import total_ordering
import re
import pandas as pd

def parse_group_string(group_string,last_mouse):
    """
    Parsing the cell meta data record from experiment.
    Input:
        group_string: string to be parsed
        last_mouse: largest mouse number in past experiments 
    """
    
    # ---RE's---
    # in l1 entry
    re_cell_prefix = re.compile(r'^\d{8}')
    re_last_cell = re.compile(r'-(\d+):')
    # in l2 entry
    re_cell_range = re.compile(r'(\d+)-(\d+)\s+')
    re_mouse_number = re.compile(r'\bm(\d+)\b')
    re_clock = re.compile(r'\bo(\d+)\b')
    re_partition = re.compile(r', ([a-z]+[0-9])\b')
    re_cell_type = re.compile(r'\bc([c,to,or,\d]+)\b')
    
    # ---Parsing---
    current_cell_prefix = None
    current_digits = None
    groups = []
    for line in group_string.split("\n"):
        if line.endswith(":"):
            # is l1 entry
            current_cell_prefix = re_cell_prefix.search(line).group(0)
            # length of last cell number
            current_digits = len(re_last_cell.search(line).group(1))
            continue
        else:
            # is L2 entry
            if current_cell_prefix == None or current_digits == None:
                print(line)
                raise ValueError("L2 entry without L1 entry.")
            try:
                start = int(re_cell_range.search(line).group(1))
                end = int(re_cell_range.search(line).group(2)) + 1
                partition = re_partition.search(line).group(1)
                cell_type = re_cell_type.search(line).group(1)
                mouse_number = int(re_mouse_number.search(line).group(1))
                clock = re_clock.search(line).group(1)
                for i in drange(start, end, current_digits):
                    # sample_name, group, partition, cell_type
                    i_sample_name = current_cell_prefix + i
                    i_group = "m" + str(last_mouse + mouse_number) + "_c" + cell_type + "_o" + clock
                    i_partition = partition
                    i_cell_type = "c" + cell_type
                    groups.append((i_sample_name, i_group, i_partition, i_cell_type))
            except AttributeError:
                print("parse_group_string: Parsing Failed - ",line)
                continue
    groups = pd.DataFrame(groups, columns = ["sample_name","group","partition","cell_type"])
    groups = groups.set_index("sample_name")
    return groups

# --- meta data: group: parser ---
norm_clock={"m22_c2_o20302100": "m22_c2_o20302100",
"m12_c2_o2130": "m12_c2_o2130",
"m13_c2_o23": "m13_c2_o2300",
"m23_c2_o22002230": "m23_c2_o22002230",
"m9_c2_o23": "m9_c2_o2300",
"m27_c2_o700730": "m27_c2_o07000730",
"m10_c2_o1": "m10_c2_o0100",
"m20_c2_o530": "m20_c2_o0530",
"m4_c2_o1112": "m4_c2_o11001200",
"m18_c2_o500530": "m18_c2_o05000530",
"m29_c2_o730": "m29_c2_o0730",
"m21_c2_o530": "m21_c2_o0530",
"m19_c2_o630700": "m19_c2_o06300700",
"m1r_c2_o11": "m1r_c2_o1100",
"m24_c2_o10": "m24_c2_o1000",
"m28_c2_o930": "m28_c2_o0930",
"m30_c2_o1330": "m30_c2_o1330",
"m31_c2_o1011": "m31_c2_o10001100",
"m26_c2_o12": "m26_c2_o1200",
"m7_c2_o10": "m7_c2_o1000",
"m1l_c2_o10": "m1l_c2_o1000",
"m3_c2_o1011": "m3_c2_o10001100",
"m2_c2_o16": "m2_c2_o1600"}

@total_ordering
class ExpTime:
    @staticmethod
    def d2(time_str):
        if len(time_str) != 2:
            raise ValueError("d2: length must be 2")
        if time_str == "00":
            return 0
        return int(time_str.lstrip("0"))
    def __init__(self,time_str):
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
            raise ValueError("time_str must be either 8 or 4 chrs")
    def __eq__(self, other):
        return (self.tah == other.tah) and (self.tam == other.tam) \
            and (self.tbh == other.tbh) and (self.tbm == other.tbm)
    def __lt__(self, other):
        #return ((self.tah < other.tah) and (self.tam < other.tam)) \
        #    or ((self.tbh < other.tbh) and (self.tbm < other.tbm))
        #return self.tah < other.tah and self.tam < other.tam
        return (self.tah*60 + self.tam) < (other.tah*60 + other.tam)
    def __str__(self):
        return str(self.tah)+":"+str(self.tam)+" "+str(self.tbh)+":"+str(self.tbm)
def add_order(annote,order):
    order = pd.Series(list(range(len(order))),index=order,name="order_index")
    return pd.concat([annote,order],axis=1)
def add_time_group_order(annote):
    # must have exp_time col
    exp_time_order = sorted(annote["exp_time"].unique(),key=ExpTime)
    exp_time_order = exp_time_order[14:] + exp_time_order[:14]
    time_group_order = pd.Series(
        list(range(len(exp_time_order))),
        index=exp_time_order,
        name="time_group_order"
    )
    return pd.merge(annote,time_group_order,left_on="exp_time",right_index=True)
def plot_cdps_time_group(annote, cdps):
    """
    Input
        cdps: contact decay profiles
        annote: sample_name as index; must have exp_time, order_index, time_group_order cols
    """
    exp_time_order = annote.sort_values("time_group_order")["exp_time"].unique()
    fig = make_subplots(
        rows=len(exp_time_order)+1,
        cols=1,
        shared_xaxes=True,
        # two types of fig have same height
        row_heights=[1/len(exp_time_order) for i in exp_time_order] + [1]
    )
    # unique return in order of appearence
    for i,t in enumerate(exp_time_order):
        subdat = annote.query('exp_time == @t')
        fig.add_trace(
            go.Scatter(
                x = subdat["order_index"],
                y = [0 for i in subdat["order_index"]],
                mode = "markers",
                name = str(ExpTime(t))
            ),
            row = i+1,
            col = 1
        )
    fig.add_trace(
        go.Heatmap(
            z = cdps.loc[annote.sort_values("order_index").index].T,
            y = cdps.columns,
            colorscale="bluyl",
            showscale=False
        ),
        row = len(exp_time_order) + 1,
        col = 1
    )
    #fig.add_trace(go.Scatter(y=[0,0,0,0],x=[1,2,3,5],mode="markers",marker_color="brown",name="o11"),row=1,col=1)
    #fig.add_trace(go.Scatter(y=[0,0,0,0],x=[2,4,5,8],mode="markers",marker_color="brown",name="o12"),row=2,col=1)
    fig.update_layout(
    {
    "plot_bgcolor": "rgba(0, 0, 0, 0)",
    "paper_bgcolor": "rgba(0, 0, 0, 0)",
    }
    )
    fig.update_xaxes(
        visible=False
    )
    fig.update_yaxes(
        visible=False
    )
    fig.update_layout(
        height = 800,
        width = 800
    )
    return fig
