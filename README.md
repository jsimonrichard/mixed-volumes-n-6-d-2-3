# Inequalities of Mixed Volumes on Six Planar Bodies in Two and Three Dimensions

This is the official repository for [Mixed volumes of zonoids and the absolute value of the Grassmannian](https://arxiv.org/abs/2404.02842), which contains all of the SageMath code used to produce the computation results in the paper.

BibTex:
```bib
coming soon...
```

## Getting Started

### Dependencies


### Two Dimensions


### Three Dimensions
```python

```


## Structure

In addition to the notebooks mentioned above, this project includes [mvol_map_problem.sage](./6_zonoids_dim_2/mvol_map_problem.sage). This file contains the `MVolMapProblem` class, which handles the parameters for the problem (number of bodies, number of dimensions, etc.) and provides cached (memoized) functions that generate important SageMath objects (like the polyhedral cone $C_{n,d}$ which is the conic hull of the image of $\Phi$). However, note that this class is unlikely to work for other dimensions.


## Contributing
Pull requests are welcome!

## License
This code has been released under the [MIT License](./LICENSE).

## Contact
If you have questions or comments feel free to email Ivan Soprunov (i.soprunov@csuohio.edu). For more code-specific questions consider reaching out to Simon Richard (j.s.richard@vikes.csuohio.edu) or Gennadiy Averkov (averkov@b-tu.de).


## Tasks
- [ ] Check for trivial proofs given duplicate triangles (sum like terms in orbit representations)
