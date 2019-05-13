# Final Project Diary
## CS 321 Spring 2019
## R. Teal Witter

### Hour 1
I looked for datasets at GEO.
I ended up finding an article with hierarchical clustering
and a dataset!!

### Hour 2
I explored the supplementary dataset to the 
article If ound.
I'm currently struggling with what data
they used to make their hierarchical clustering
figure.

### Hour 3
I figured out what data they used for Figure 2 b :)
It's the first tab in the first supplementary excel file.

### Hour 4
I worked on reading the data into my python file.
I converted to a CSV then searched for the genes used in Fig 2b.

### Hour 5
Extracted distance from data!
Ran hierarchical clustering on the data and close-ish results :)

### Hour 6
I gained familiarity with the SciPy dendrogram function.
I looked at your code and the example online

### Hour 7
I refactored my code to create a linkage that works as input
to the dendrogram function.
I made several dendrograms!

### Hour 8
I figured out that I need to cluster the tumor for the heat map 
as well.
I extracted their IDs and plotted them with the seaborn clustermap
but now I need to figure out how to combine a heat map and clustering.

### Hour 9
I researched different possible clustering techniques to replicate
Figure 2b: none of single, complete, average, weighted, centroid,
median, or ward exactly replicate the clustering.
The paper talks about ``unsupervised'' clustering so that might
be what they used.

### Hour 10
I did a `spike' into what it would take to create a clustermap.
Making it on my own seems way too hard so I'm going to split it
up into clustering of the genes, the tumors and then a heatmap.

### Hour 11
I started comparing the tumor clustering with their PAM50 identifications.
The complete and average are remarkably accurate but it took awhile
to find the right threshold for each one: complete is about 10 while
average is about 8.

### Hour 12
I began working on the poster
I did all the initial formatting and layout issues and started
brainstorming what figures would be most illustrative with
the limited poster size.

### Hour 13
I worked on generating accurate graphs for the clusterings.
Couldn't figure out a CS way to copy the flat clusters over so I 
manually did it RIP :(

### Hour 14
Realized I didn't tranpose before finding clustering of supposed
"tumor" clustering. I had to manually repeat the process of
copying over the flat clusters.

### Hour 15
Finished the algorithm section of my poster.
Tried to color PAM50 analysis, can't find a way to automate it.

### Hour 16
I decided to print out simplified versions of the PAM50 subtypes
and manually color.
I refactored and commented code. (Yay!)

### Hour 17
Wrote goal, outline for background and experiment.
Worked on layout of poster (you're right, this is insanely time
comsuming).

### Hour 18
Worked on exporting heatmap. Lined up the colors to be the same!

### Hour 19
Put the clustered heatmaps together.
I have a new appreciation for abstraction: it took me 2 cluster images,
one (manually colored) string of PAM50 subtypes, 
and two heatmaps to make each one of the three `clustermaps'.

### Hour 20
Put all the clustermaps on the powerpoint.
Wrote out all the words and arranged things!

### Hour 21
I finalized the poster and showed it to friends.
Practiced introducing the topics and presenting.

### Hour 22
I made the PAM50 subtypes more visible.
I also reorganized the poster to be less cluttered.
