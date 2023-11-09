"""run with 'plotting' conda env """
#  import statements
import pandas as pd
from matplotlib_venn import venn2
from matplotlib import pyplot as plt

#  create venn diagrams

## https://www.geeksforgeeks.org/how-to-create-and-customize-venn-diagrams-in-python/

## venn #1
venn2(subsets=(69, 127, 17), set_labels=('melanoma', 'basal cell carcinoma'))
#plt.suptitle("Tumor-Tumor Overlap -- Upregulated TFs in Macrophage", fontsize=18)  # https://stackoverflow.com/questions/1388450/giving-graphs-a-subtitle-in-matplotlib
plt.title("TFs Up-Regulated in Macrophage \n adjusted p-value: 0.041")
plt.show()

## venn #2
venn2(subsets=(133, 78, 26), set_labels=('melanoma', 'basal cell carcinoma'))
plt.title("TFs Down-Regulated in Macrophage \n adjusted p-value: 0.00034")
plt.show()

## venn #3
venn2(subsets=(187, 47, 52), set_labels=('melanoma', 'basal cell carcinoma'))
plt.title("TFs Up-Regulated in CD8 T Cells \n adjusted p-value: 3.65E-16")
plt.show()

## venn #4
venn2(subsets=(103, 212, 63), set_labels=('melanoma', 'basal cell carcinoma'))
plt.title("TFs Down-Regulated in CD8 T Cells \n adjusted p-value: 3.62E-8")
plt.show()

## venn #5
venn2(subsets=(31, 48, 6), set_labels=('melanoma', 'basal cell carcinoma'))
plt.title("TFs Up-Regulated in NK Cells \n adjusted p-value: 0.019")
plt.show()

## venn #6
venn2(subsets=(49, 90, 25), set_labels=('melanoma', 'basal cell carcinoma'))
plt.title("TFs Down-Regulated in NK Cells \n adjusted p-value: 7.13E-10")
plt.show()

