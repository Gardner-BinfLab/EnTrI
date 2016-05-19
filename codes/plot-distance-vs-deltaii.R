require(ggplot2)
require(stringr)
# require(dplyr)
diipath <- "../results/duplication_distance_calculator/deltaiis/"
dspath <- "../results/duplication_distance_calculator/distances/"
list_of_files <- list.files(path=dspath, pattern="*.txt", full.names=FALSE, recursive=FALSE)
dii <- c()
ds <- c()
genes <- c()
for (filename in list_of_files)
{
  dii_tbl <- read.table(paste(diipath,'/',filename,sep = ''), sep = '\t', header = FALSE)
  ds_tbl <- read.table(paste(dspath,'/',filename,sep = ''), sep = '\t', header = FALSE)
  if (nrow(dii_tbl) > 1)
  {
    for (i in seq(1,(nrow(dii_tbl)-1)))
    {
      for (j in seq((i+2), ncol(dii_tbl)))
      {
        genei = str_match(dii_tbl[i,1], "([[:graph:]]*)\\_[[:digit:]]+\\-[[:digit:]]+")[2]
        genomei = str_match(genei, "([[:alnum:]]*)\\_[[:graph:]]+")[2]
        if (is.na(genomei))
        {
          genomei = str_match(genei, "([[:alpha:]]*)[[:digit:]]+")[2]
        }
        genej = str_match(dii_tbl[j-1,1], "([[:graph:]]*)\\_[[:digit:]]+\\-[[:digit:]]+")[2]
        genomej = str_match(genej, "([[:alnum:]]*)\\_[[:graph:]]+")[2]
        if (is.na(genomej))
        {
          genomej = str_match(genej, "([[:alpha:]]*)[[:digit:]]+")[2]
        }
        if (dii_tbl[i,j] != -1 & genomei  == genomej)
        {
          dii <- c(dii, dii_tbl[i,j])
          ds <- c(ds, ds_tbl[i,j])
          genes <- rbind(genes, c(genei, genej))
        }
      }
    }
  }
}

pdf("../results/distance-vs-deltaii.pdf")

# plot(ds, dii, ylim=c(0,2))
# lines(lowess(ds, dii, f = .2), col = 2, lwd=5)
# regln <- lm(dii~ds)
# abline(regln, col="red")

ds_quartiles <- c()
thresholds <- quantile(ds)
for (value in ds)
{
  if (value < thresholds[2])
  {
    ds_quartiles <- c(ds_quartiles, 'New')
  }
  else if (value < thresholds[3])
  {
    ds_quartiles <- c(ds_quartiles, 'Older new')
  }
  else if (value < thresholds[4])
  {
    ds_quartiles <- c(ds_quartiles, 'Newer old')
  }
  else
  {
    ds_quartiles <- c(ds_quartiles, 'Old')
  }
}
df <- data.frame(ds_quartiles, dii)
df$ds_quartiles <- factor(df$ds_quartiles,levels = c("New", "Older new", "Newer old", "Old"),ordered = TRUE)
ggplot(df)+geom_violin(aes(x=ds_quartiles,y=dii,group=ds_quartiles),color="black",fill="grey",size=2)+
  theme(axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"))+   # build the base violins
  labs(x="Distance", y = "Delta ii")
dev.off()
#http://stackoverflow.com/questions/22278951/combining-violin-plot-with-box-plot