
#clear the R environment
rm(list = ls())
# setting working directory
setwd("C:/Users/Babatunde_Oriyomi1/Documents/Babs_phylo")

#install new packages
install.packages("ape")
install.packages("phangorn")
install.packages("seqinr")
install.packages("ggtree")
Install.packages("tidyverse")

# import libraries
library(ape)
library(seqinr)
library(phangorn)
library(tidyverse)
library(ggtree)


# URL of the raw FASTA file on GitHub
fasta_url <- "https://raw.githubusercontent.com/Teespice/Babatunde_Tutorial-/main/sample_data.fasta"

# Local file path where you want to save the downloaded file
local_file_path <- "sample_data.fasta"

# Download the file
download.file(url = fasta_url, destfile = local_file_path, method = "auto")

# Check if the file has been downloaded successfully
if (file.exists(local_file_path)) {
  cat("FASTA file downloaded successfully.")
} else {
  cat("Failed to download the FASTA file.")
}

# Read the downloaded file as text
fasta_text <- readLines("sample_data.fasta")

# Print the contents to the console
cat(fasta_text, sep = "\n")

# Read the FASTA file
fasta_file <- "sample_data.fasta"
#view fasta file
file.show("sample_data.fasta")

# reading and aligning sequences
DNAseqs <- clustal(read.dna("sample_data.fasta", format = "fasta"))

# converting to phydat
sequence_phydat <- phyDat(DNAseqs, type = "DNA", levels = NULL)

# testing models
mt <- modelTest(sequence_phydat)
mt <- mt[order(mt$BIC),]

# Compute the distance matrix
dna_dist <- dist.ml(sequence_phydat, model="JC69")

# Generate the phylogenetic tree using the UPGMA method
mytree1 <- NJ(dna_dist)

plot(mytree1)    # Distance matrix tree, but not rooted at midpoint.

# saving an image of the distance matrix tree that is rooted at the midpoint
png(filename = "DISTMATRIXTREE.png",width = 10, height = 10, units = "in", res = 300)
plot(midpoint(mytree1), main = "Distance Matrix tree")
title("Phylogenetic tree for mosquito species", line = -44)
add.scale.bar()
dev.off()

# Save unrooted tree
png(filename = "unrootedTREE.png",width = 10, height = 10, units = "in", res = 300)
plot(mytree1, "unrooted")
dev.off()

# plotting ML tree
fit <- pml(mytree1, sequence_phydat)
print(fit)
fitJC <- optim.pml(fit, model = "JC")
logLik(fitJC)
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))

## This is the actual ML tree:
MLtree <- fitJC$tree
# Plot
print(MLtree)

# Root this one at the midpoint also and plot again
plot(midpoint(MLtree), main = "Maximum Likelihood Tree")

## The ML tree is almost the same as the distance matrix tree, but it has A. meras in a clade with 4 others.
## So maybe this shows that the distance matrix is a pretty good fit in this case, do you think?

## Export it:
png(filename = "ML_TREE.png",width = 10, height = 10, units = "in", res = 300)
plotBS(midpoint(MLtree), main = "Maximum Likelihood Tree", bs)
title("Phylogenetic tree for mosquito species", line = -44)
dev.off()


## You can also write these to FigTree files if you want to work with them further.
# Write to FigTree file
write.tree(midpoint(MLtree), file="bootstrap_example.tre")


## Read tree into ggtree for aesthetics:

#since the tree has negative edge lengths and ggtree will complain about it.
options(ignore.negative.edge=TRUE)

## Plot distance matrix tree:
q <- ggtree(midpoint(mytree1)) +
  geom_tiplab(size=2.6, align=FALSE, linesize=.5)  
print(q)


# use ggsave to save the plot.
# Sometimes you have to adjust the font size (in the geom_tiblab above), width, and height because it looks different in the image.
ggsave(filename = "test_tree.png", plot = q, width = 12, height = 8, dpi = 300, units = "in")



## Plot maximum liklihood tree:
maxtree <- ggtree(midpoint(MLtree)) +
  geom_tiplab(size=2.6, align=FALSE, linesize=.5) +
  labs(caption = "Maximum Likelihood Tree") +
  theme(plot.caption = element_text(hjust = 0.5, size = rel(.75))) 
print(maxtree)


# Also, it looks like we can view all the trees in the bootstrap list that was generated like this,
# if you ever wanted to.
b_tree <- ggtree(bs) + facet_wrap(~.id, ncol=10)
print(b_tree)


## To add colors, first we need to know the tip labels (the taxa).
print(MLtree$tip.label)

# Usually there is data in the fasta file that lets us match colors to something, but there isn't in this case.
# So we will need to match a color manually to each one of these tip labels. 

# Define your colors in variables, and you can switch to whatever colors you like.
color0 <- "grey"  # red
color1<- "blue"  # orange
color2 <- "green"  # Green
color3 <- "orange"  # Blue


# Create a data frame that manually maps the tip labels to the colors. We know the tip labels are in
# MLtree.tip.label, so we can use that here. The colors have to be manually added to the correct places.
tip_data <- data.frame(
  tip.label = MLtree$tip.label,
  color = c(color0, color0, color1, color1, color0, color2, color2, color2, color2, color2,
            color3, color2, color2, color2, color3, color3, color3, color3, color3)
)


# First make a basic plot of the tree without colors, since the colors get added in the step below
p <- ggtree(midpoint(MLtree)) + 
  geom_tiplab(size = 2.6, align = FALSE)
print(p)


# Add the tip data we created to the tree; the "color" referenced below is the value in tip_data
p <- p %<+% tip_data + 
  geom_tippoint(aes(color = color), size = 2)  

# If you plot it like this, it will give you colors, but it won't use the colors specified. (That's the next step.)
print(p)


# If we use scale_color_identity, it will give the exact colors you ask for in the color variables.
# For me, this step also removes the inapplicable color legend.
p <- p + scale_color_identity()
print(p)


# You might want a scale bar that represents substitutions per site. (I think that a difference of 0.10 on the scale means there is a 10% difference between clades.)
p <- p + geom_treescale(x = 0.1, y = 0.1, width = 0.1, fontsize = 3, offset = 0.02, color = "black")
print(p)


## All in one step:
s <- ggtree(midpoint(MLtree)) %<+% tip_data + 
  geom_tiplab(size = 2.6, align = FALSE) + title(main = "Maximum Likelihood tree")+
  geom_tippoint(aes(color = color), size = 2.5)  +   ## < Change the size of colored tips here
  geom_treescale(x = 0.1, y = 0.1, width = 0.1, fontsize = 3, offset = 0.02, color = "black") +
  scale_color_identity()
print(s)

ggsave(filename = "colored_tree.png", plot = s, width = 16.7, height = 8, dpi = 300, units = "in")


