# -*- coding: utf-8 -*-
"""
Created on Thu May 28 13:08:31 2020

@author: Yuni
"""

from skimage import data, segmentation, filters, color
from skimage.future import graph
from matplotlib import pyplot as plt


def weight_boundary(graph, src, dst, n):
    """
    Handle merging of nodes of a region boundary region adjacency graph.

    This function computes the `"weight"` and the count `"count"`
    attributes of the edge between `n` and the node formed after
    merging `src` and `dst`.


    Parameters
    ----------
    graph : RAG
        The graph under consideration.
    src, dst : int
        The vertices in `graph` to be merged.
    n : int
        A neighbor of `src` or `dst` or both.

    Returns
    -------
    data : dict
        A dictionary with the "weight" and "count" attributes to be
        assigned for the merged node.

    """
    default = {'weight': 0.0, 'count': 0}

    count_src = graph[src].get(n, default)['count']
    count_dst = graph[dst].get(n, default)['count']

    weight_src = graph[src].get(n, default)['weight']
    weight_dst = graph[dst].get(n, default)['weight']

    count = count_src + count_dst
    return {
        'count': count,
        'weight': (count_src * weight_src + count_dst * weight_dst)/count
    }


def merge_boundary(graph, src, dst):
    """Call back called before merging 2 nodes.

    In this case we don't need to do any computation here.
    """
    pass


name='1'+"_BF.tif"  
img = cv2.imread(name)
edges = filters.sobel(color.rgb2gray(img))
labels = segmentation.slic(img, compactness=30, n_segments=800, start_label=1)
g = graph.rag_boundary(labels, edges)

graph.show_rag(labels, g, img)
plt.title('Initial RAG')

labels2 = graph.merge_hierarchical(labels, g, thresh=0.0001, rag_copy=False,
                                   in_place_merge=True,
                                   merge_func=merge_boundary,
                                   weight_func=weight_boundary)

graph.show_rag(labels, g, img)
plt.title('RAG after hierarchical merging')

plt.figure()
out = color.label2rgb(labels2, img, kind='avg', bg_label=0)
plt.imshow(out)
plt.title('Final segmentation')

plt.show()











import matplotlib.pyplot as plt

from skimage.filters import sobel
from skimage.measure import label
from skimage.segmentation import slic, join_segmentations, watershed
from skimage.color import label2rgb
from skimage import data

coins = cv2.imread("1_BF.tif_nucle_mask.jpg")
coins = cv2.cvtColor(coins, cv2.COLOR_BGR2GRAY)
plt.imshow(coins)
# Make segmentation using edge-detection and watershed.
edges = sobel(coins)
plt.imshow(edges)
# Identify some background and foreground pixels from the intensity values.
# These pixels are used as seeds for watershed.
markers = coins
foreground, background = 1, 2
markers[coins < 1] = background
markers[coins > 100.0] = foreground
plt.imshow(coins < 30.0)
plt.imshow(foreground)
ws = watershed(edges, markers)
seg1 = label(ws == foreground)
plt.imshow(seg1)

# Combine the two.
segj = join_segmentations(seg1, seg2)




