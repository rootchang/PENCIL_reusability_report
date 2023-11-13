# PENCIL\_reusability\_report
This is a reusability report of [PENCIL](https://doi.org/10.1038/s42256-023-00656-y).

<p align="center">
  <img src="./img/schema.jpg" width = "1000" alt="method" align=center />
</p>

## 1. How to reproduce figures in the reusability report
To reproduce the figures in the reusability report, please run the R markdown scripts (.Rmd) in this repository.

## 2. How to train PENCIL on training data and apply it on new data

I have made a pipeline `Run_PENCIL_csv.py` for running PENCIL with the .csv input files. The basic usage is as follows:

```{bash}
python Run_PENCIL_csv.py p1 p2 p3 p4 p5 p6 p7 p8
```

Table 1. Explanation for the parameters used by `Run_PENCIL_csv.py`.

| Parameter | Description | Default |
|----------|----------|----------|
| p1   | Input scaled gene expression .csv file | required, no default |
| p2   | Input cell phenotypic label .csv file | required, no default |
| p3   | embedding .csv file | required, no default |
| p4   | data name, folder in the results directory to store result | required, no default |
| p5   | experiment id, PENCIL will make a new folder in the results directory use this name | required, no default |
| p6   | PENCIL mode, can be either 'multi-classification' or 'regression' | required, no default |
| p7   | class_weights, separated by ',' | required, no default |
| p8   | Input scaled gene expression .csv files (full path and file name) for test data, separated by ',' | required, no default |

## Citation


## Contact
changtiangen@gmail.com

