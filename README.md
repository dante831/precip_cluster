# precip_cluster

A code suite to calculate precipitation clusters and their distribution on a power-law plot. 
Running this code together with data downloaded from this [link](https://www.dropbox.com/sh/c4faexb4dax40xz/AABPnFfP6Tm9IIimn9T5VzwXa?dl=0) 
reproduces Fig. 2 in [O'Gorman et al, 2021](https://royalsocietypublishing.org/doi/10.1098/rsta.2019.0543). 

## How to run: 

### Download data

Download data from this link: https://www.dropbox.com/sh/c4faexb4dax40xz/AABPnFfP6Tm9IIimn9T5VzwXa?dl=0. 
You may have to download the data from three folders separately since they are too large to download at once. 
Make a new folder in the root directory named "data". Then, put the downloaded data in its original structure
in the folder. 

### Perform analysis

In the root directory, open a terminal and type the following code: 

    python main.py
    
The produced figure will be in a folder named "figures/"

## License

This code is distributed under the [MIT](http://opensource.org/licenses/mit-license.php) license
