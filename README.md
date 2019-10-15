## Needle
Needle provides a space-efficient data structure to index a large amount of NGS data and allows fast searches through these indices.
Due to the space-efficiency of one index, it is affordable to create multiple indices with different expression rates. Therefore, a semi-quantitative analysis of the data becomes possible. Needle is based on Interleaved Bloom Filters, which is a compact and efficient structure to store multiple Bloom Filters. Furthermore, Needle uses a windowing scheme (also called Minimizers) to reduce the amount of data to store.  

## Build

Needle is depending on the seqan3 library, at the moment it is necessary to use Enrico Seiler's branch "" (), where the IBF is implemented. Soon, this branch should be included in the seqan3 library.
Assuming seqan3 can be found in "${CMAKE_SOURCE_DIR}/../", Needle can be build following these commands:

```
git clone 
mkdir build-needle & cd build-needle
cmake ../needle
make
```

## Create an IBF


```
needle-ibf
```

## Search


```
needle-search
```
