# extract code cells from jupyter ipynb, use when parsing fails
# Usage:
#   python rescue_jupyter.py output input
# Change C if you want to extract more code cells
import re
import sys
output, input_ = sys.argv[1], sys.argv[2]
def find_code(text):
    pattern = r'\"cell_type\":\s?\"code\".*\n(.*\n)*?\s*\"source\":\s*\[((.*\n)*?)\s*\]\n'
    # (.*\n) means any new line
    # python code signature:
    # starts with "source": [
    #    '\s*\"source\":\s*\['
    # ends with single ] line
    #    '\s*\]\n'
    # (.*\n)*? make 2  "any new line" pattern non-greedy, critical, otherwise it will match to the last line
    res = re.findall(pattern, text)
    if res:
        return (i[1] for i in res)
def print_code(source):
    # strip \n and " from each line
    noblank = (i.rstrip(",") for i in (i.strip() for i in source.split("\n")))
    body = (i.rstrip("\\n") for i in  (i.strip("\"") for i in noblank))
    return body
with open(input_, "r") as f:
    with open(output, "w") as f2:
        C = 1000 # max number of code cells to extract, use for debugging
        for i, code in enumerate(find_code(f.read())):
            for line in print_code(code):
                f2.write(line+"\n")
            if i > C:
                break