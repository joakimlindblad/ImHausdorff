# Hausdorff
Matlab code for computing Hausdorff distance between binary images or segmentation maps.

---

Noticing that all Hausdorff distance functions I could find on MATLAB File Exchange
seem to compute all pairwise point-to-point distances making them very slow for large
images, I felt compelled to write a version which uses bwdist and thus scale linearly
for large image data. The code should work for images of any dimension.
Do try the examples in the help text.

Output of 'euclidean_precise' is equal to the below slower point-cloud version

```matlab
if isempty(a) || isempty(b)
        if isempty(a) && isempty(b)
                HD=0;
        else
                HD=inf;
        end
        return
end

% Max of dist from A to B and B to A
%D = pdist2(a,b,'squaredeuclidean'); 
%HD = sqrt(max(max(min(D,[],1)),max(min(D,[],2)))); 

% Less memory hungry version
D1 = pdist2(a,b,'squaredeuclidean','Smallest',1);
D2 = pdist2(b,a,'squaredeuclidean','Smallest',1);
HD = sqrt(max(max(D1),max(D2)));
```

/ Joakim
