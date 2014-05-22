WBHC-Project
============

Code for Word-Based Heirarchical Classification model

This project develops a simple algorithm to make accurate feature-based class predictions when your classes are highly unbalanced AND they are heirarchically structured. It works by recursively splitting the dataset such that: a) the class heirarchy is preserved, and b) the split produces two groups which are as similar in size as possible, and c) the groups contain data which is 'similar' in its features (so the most difficult classifications occur last). Even though using binary classification in a mult-class problem is sub-optimal, the gains in data balance more than make up for it (data balance is one of the most difficult issues for machine learning currently). 

Once these binary splits are determined, each binary classification is performed using a randomForest algorithm (with regularization), which can be done in parallel. For classifying unknown data, the data are runs through all models, and then the prediction probabilities from each model are chained together using a tree-like structure, eventually resulting in a prediction probability for assigning an unknown sample to each taxonomic group along the entire heirarchy. 

I use this model to classify genetic sequences into heirarchically structured taxonomic classes (e.g. Order, Family, Genus), using word (or K-mer) 'waiting times' as a vector of features for each sequence. More specifically, I use the mean and variance of the gaps between the reoccurence of a word (the waiting times in Poisson parlance), for all possible words of size 1 to k to form a large feature vector describing the statistical properties of that sequence (vector length = sum<sub>1:k</sub>(2*(4<sup>k</sup>))). 
