__author__ = 'sam.royston'
import numpy as np
import matplotlib.pyplot as plt
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

class PartialTree(object):
    def __init__(self, height, sample_density=0.01, max_samples = 10000):
        # start with a single branch
        self.branch = np.random.normal(mu,sigma,height+1)
        self.branch[0] = 0
        self.max = {"sum":0, "path":self.branch}
        self.min = {"sum":0, "path":self.branch}
        self.minima = []
        self.maxima = []
        self.max_samples = max_samples
        for i in xrange (1, height + 2):
            self.minima.append({"sum":0, "path":np.zeros(i)})
            self.maxima.append({"sum":0, "path":np.zeros(i)})
        self.sampling = []
        self.sample_density = sample_density

    def check_minmax(self, density = 1):
        """
        check current branch against past min/max
        """
        i = 1
        for minimum in self.minima[::density]:
            range = self.branch[0:i]
            path_sum = np.sum(range)
            if path_sum < minimum["sum"]:
                self.minima[i-1] = {"sum":path_sum, "path":np.copy(self.branch[0:i])}
            i += density
        i = 1
        for maximum in self.maxima[::density]:
            range = self.branch[0:i]
            path_sum = np.sum(range)
            if path_sum > maximum["sum"]:
                self.maxima[i-1] = {"sum":path_sum, "path":np.copy(self.branch[0:i])}
            i += density

    def pick_at_random(self):
        """
        update random selection of branches for viz if needed
        """
        if random() < self.sample_density and len(self.sampling) < self.max_samples:
            self.sampling.append(np.copy(self.branch))

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

    ## modify the branch for each path enumerated using the binary representation above
    while int(string_state_representation,2) >= 0:
        new_state_representation = bin(int(string_state_representation,2) - 1)

        ## this trick tells us how far up the branch to replace
        string_indices_to_change = bin(int(new_state_representation,2) ^ int(string_state_representation,2))
        string_state_representation = new_state_representation

        ## how long is the XOR ? ... (- 2) because of '0b' chars in string representation
        deepness = len(string_indices_to_change) - 2
        index = len(partial_tree.branch) - deepness
        partial_tree.branch[index:] = np.random.normal(mu,sigma,deepness)

        ## checks for new extremes ~ costly b/c default param looks for extremes at every level of the tree
        partial_tree.check_minmax()

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

def take_average(n, depth):
    """
    take average values in n independent GREM trials
    """
    averages_min = [{"path":np.zeros(i), "sum":0} for i in xrange(1,depth + 2)]
    averages_max = [{"path":np.zeros(i), "sum":0} for i in xrange(1,depth + 2)]
    for t in xrange(0,n):
        c = '\033[92m' if t == n - 1 else '\033[91m'
        sys.stdout.write("\r" + c + str((t * 1.0/n * 1.0) * 100) + '% \033[0m' )
        sys.stdout.flush()
        results = permute_branches(DEPTH)
        for i,m in enumerate(results.minima):
            averages_min[i]["path"] = averages_min[i]["path"] + (m["path"] / n)
        for i,m in enumerate(results.maxima):
            averages_max[i]["path"] = averages_max[i]["path"] + (m["path"] / n)

    return averages_min, averages_max

## take averages
_averages_min, _averages_max = take_average(averaging_time, DEPTH)

## draw to matplotlib environment
"""
b: blue
g: green
r: red
c: cyan
m: magenta
y: yellow
k: black
w: white

Gray shades can be given as a string encoding a float in the 0-1 range, e.g.:

color = '0.75'
For a greater range of colors, you have two options. You can specify the color using an html hex string, as in:

color = '#eeefff'
"""

colors = "bgrcm"
col =0
for minimum,maximum in zip(_averages_min, _averages_max):
    col += 1
    draw_array(minimum["path"], colors[col % len(colors)] + "-", alpha=0.2)
    draw_array(maximum["path"], colors[col % len(colors)] + "-", alpha=0.2)
plt.axis('off')
plt.savefig("../GREM.pdf")
plt.show()







