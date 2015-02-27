# Median scale data

# Notes:
# 1. 

# Input:
# data (positive values: raw intensities or pseudo-counts)

# Output:
# Median scaled data

medScal = function (data) {
  meds = apply(data, 2, median, na.rm=TRUE)
  sweep(data, 2, meds, "/")
}