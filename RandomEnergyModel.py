__author__ = 'sam.royston'
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys

from random import random

DEPTH = int(sys.argv[1])
if len(sys.argv) > 2:
    averaging_time = int(sys.argv[2])
else:
    averaging_time = 1

canvas_height = 800
canvas_width = 1400
scaling_factor = 1
start_pos = 0
mu,sigma = 0, 1
pp = PdfPages('GREM.pdf')

class PartialTree(object):
    def __init__(self, height, sample_density=0.01, max_samples = 10000):
        # start with a single branch
        self.branch = np.random.normal(mu,sigma,height+1)
        self.branch[0] = 0
        self.max = {"sum":-1000, "path":self.branch}
        self.min = {"sum":1000, "path":self.branch}
        self.max_norm_l1 = {"norm":0, "path":self.branch}
        self.min_norm_l1 = {"norm":1000000, "path":self.branch}
        self.minima = []
        self.maxima = []
        self.max_samples = max_samples
        for i in xrange (1, DEPTH + 1):
            self.minima.append({"sum":1000, "path":self.branch})
            self.maxima.append({"sum":-1000, "path":self.branch})
        self.sampling = []
        self.sample_density = sample_density

    def check_minmax(self):
        """
        check current branch against past min/max
        """
        i = 1
        for minimum in self.minima:
            range = self.branch[0:i]
            path_sum = np.sum(range)
            if path_sum < minimum["sum"]:
                self.minima[i - 1] = {"sum":path_sum, "path":np.copy(self.branch[0:i])}
            i += 1
        i = 1
        for maximum in self.maxima:
            range = self.branch[0:i]
            path_sum = np.sum(range)
            if path_sum > maximum["sum"]:
                self.maxima[i - 1] = {"sum":path_sum, "path":np.copy(self.branch[0:i])}
            i += 1

    def replace_path_at_height(self,h):
        """
        replace all edges below height - h
        """
        self.branch[self.branch.size - h:] = np.random.normal(mu,sigma,h)

    def fuse(self,tree):
        """
        merge with other tree
        """
        self.branch = np.hstack( zip(self.branch,tree.branch) )

    def sort_samples(self, h):
        """
        sort samples by their sum
        """
        self.sampling = sorted(self.sampling, key = lambda x:x["sum"])

    def permute_below(self, h):
        """
        recursively try all branches below a certain level
        """
        sum = np.sum(self.branch)
        if sum > self.max["sum"]:
            self.max = {"sum":sum, "path":self.branch}
        elif sum <= self.min["sum"]:
            self.min = {"sum":sum, "path": self.branch}
        return self.permute_below(h - 1)

    def print_structure(self):
        """
        show hierarchy
        """
        pass

def permute_branches(height):
    """
    simulate all possible paths down our tree
    """
    string_state_representation = bin(2**height - 1)
    partial_tree = PartialTree(height)
    while int(string_state_representation,2) >= 0:
        new_state_representation = bin(int(string_state_representation,2) - 1)
        string_indices_to_change = bin(int(new_state_representation,2) ^ int(string_state_representation,2))
        string_state_representation = new_state_representation
        deepness = len(string_indices_to_change) - 2
        index = len(partial_tree.branch) - deepness
        partial_tree.branch[index:] = np.random.normal(mu,sigma,deepness)
        path_sum = np.sum(partial_tree.branch)
        norm_l1 = np.linalg.norm(partial_tree.branch, ord=1)

        if random() < partial_tree.sample_density and len(partial_tree.sampling) < partial_tree.max_samples:
            partial_tree.sampling.append(np.copy(partial_tree.branch))
        partial_tree.check_minmax()
        if path_sum > partial_tree.max["sum"]:
            partial_tree.max = {"sum":path_sum, "path":np.copy(partial_tree.branch)}
        elif path_sum <= partial_tree.min["sum"]:
            partial_tree.min = {"sum":path_sum, "path": np.copy(partial_tree.branch)}
        if norm_l1 < partial_tree.min_norm_l1["norm"]:
            partial_tree.min_norm_l1 = {"norm":norm_l1, "path": np.copy(partial_tree.branch)}
        if norm_l1 > partial_tree.max_norm_l1["norm"]:
            partial_tree.max_norm = {"norm":norm_l1, "path": np.copy(partial_tree.branch)}

    return partial_tree



def draw_array(a, color, width = 1, alpha = 0.005):
    """
    draw an array of energies
    """
    l = a.tolist()[1:]
    y = [0]
    y.extend([np.sum(a[0:l.index(element)+2]) for element in l])
    x = xrange (0, len(l) + 1)
    plt.plot(x,y,color, linewidth=width, alpha=alpha)
    return x,y

colors = "bgrcmy"
col =0
results = permute_branches(DEPTH)
averages_min = results.minima
averages_max = results.maxima

for t in xrange(0,averaging_time):
    print str((t * 1.0/averaging_time * 1.0) * 100) + "%"
    results = permute_branches(DEPTH)
    for i,m in enumerate(results.minima):
        averages_min[i]["path"] = averages_min[i]["path"] + m["path"]
    for m in averages_min:
        m["path"] = m["path"] 

    for i,m in enumerate(results.maxima):
        averages_max[i]["path"] = averages_max[i]["path"] + m["path"]
    for m in averages_max:
        m["path"] = m["path"]

for minimum,maximum in zip(averages_min, averages_max):
    col += 1
    col %= 6
    draw_array(minimum["path"], colors[col] + "-", alpha=0.2)
    draw_array(maximum["path"], colors[col] +  "-", alpha=0.2)

# for sample in results.sampling:
#     draw_array(sample,"k-",width=1)
# plt.plot([0, DEPTH], [0, np.sum(results.max["path"])],"k-", linewidth=3.0, linestyle='dashed')
# plt.plot([0, DEPTH], [0, np.sum(results.min["path"])],"k-", linewidth=3.0, linestyle='dashed')
print "________"



#
# X,min_y = draw_array(results.min["path"], "k-", alpha=0.3)
# X,max_y = draw_array(results.max["path"], "k-", alpha=0.3)
# draw_array(results.min_norm_l1["path"], "r-", alpha=0.8)
# draw_array(results.max_norm_l1["path"], "b-", alpha=0.8)
# new_line_min = np.interp(X, [0, DEPTH], [0, np.sum(results.min["path"])])
# new_line_max = np.interp(X, [0, DEPTH], [0, np.sum(results.max["path"])])
#
# plt.fill_between(X, min_y, new_line_min, where=min_y<=new_line_min, facecolor='teal', alpha = 0.1, interpolate=True)
# plt.fill_between(X, max_y, new_line_max, where=max_y>=new_line_max, facecolor='teal', alpha = 0.1, interpolate=True)
# plt.fill_between(X, min_y, new_line_min, where=min_y>=new_line_min, facecolor='yellow', alpha = 0.1, interpolate=True)
# plt.fill_between(X, max_y, new_line_max, where=max_y<=new_line_max, facecolor='yellow', alpha = 0.1, interpolate=True)

plt.axis('off')
plt.savefig("GREM.pdf")

plt.show()







