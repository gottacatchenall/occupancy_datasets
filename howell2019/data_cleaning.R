coords = read.csv("./data/coords.dryad.csv")
colnames(coords) = c("site", "x","y")

water_level = read.csv("./data/water.dryad.csv")
colnames(water_level) = c("site", "water_level")
head(water_level)

# converting their matrix style to a dataframe i can use
mat_to_df = function(filename){
    pre = read.csv(filename)
    data_matrix = array(data.matrix(pre[,2:45], 47), dim=c(47, 3, 15))   # sites x sample attempt x year array 
    df = data.frame(matrix(ncol=4))
    colnames(df) = c("site", "year", "sample", "val")
    row = 1
    for (site in seq(1,47)){
        for (year in seq(1,15)){
            val = 0
            for (attempt in seq(1,3)){      
                df[row,] = c(site, year, attempt, data_matrix[site, attempt, year])
                row = row + 1   
            }
        }
    }
    return(df)
}

# occupancy
occupancy = mat_to_df("./data/y.wide.dryad.csv")
colnames(occupancy) = c("site", "year", "sample", "occupancy")

# wind 
wind = mat_to_df("./data/wind.wide.dryad.csv")
colnames(wind) = c("site", "year", "sample", "wind")


temp = mat_to_df("./data/temp.wide.dryad.csv")
colnames(temp) = c("site", "year", "sample", "temp")


library(tidyverse)

full_df = coords %>% 
        left_join(water_level) %>% 
        left_join(wind) %>% 
       left_join(temp)  %>%  
       left_join(occupancy) %>% filter(is.na(occupancy) == FALSE)

write.csv(full_df, file="howell2019.csv")