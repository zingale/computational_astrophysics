# Using These Notes

These notes are built via [Jupyter book](https://jupyterbook.org/), as
a collection of [Jupyter](https://jupyter.org/) notebooks and markdown
pages.

The course is on Github at:
https://github.com/zingale/computational_astrophysics, and the course
website is built automatically via a Github action each time a change
is pushed.

```{tip}
If you find any problems or have suggestions for improving the notes,
feel free to create an [issue](https://github.com/zingale/computational_astrophysics/issues)  or [pull request](https://github.com/zingale/computational_astrophysics/pulls) at the Github repo.
```

````{note}
These notes make use of [MyST markdown](https://mystmd.org/) and [sphinx{design}](https://sphinx-design.readthedocs.io/en/furo-theme/index.html) for some features.  If you download and run the notebooks
on your own, for the best experience you should install these packages, which can be done as:

```
pip install jupyterlab_myst sphinx_design
```
````

## Interactive Usage

For the Jupyter notebooks in this collection, there are a few ways to
access them to run them on your own.

* clicking on the {octicon}`download` icon in the upper right let's
  you download the raw notebook so you can run it on your local
  computer.

* clicking on the {octicon}`rocket` icon in the upper right will allow
  you to run the notebook directly in the cloud.  There are 2 different
  compute clouds:

  * [mybinder](https://mybinder.org/) : this is an open project with
    ties to the Jupyter project.  It can take a few minutes for the
    page to appear if it hasn't been accessed recently, but then it
    will give you the standard Jupyter experience.

  * [Google colab](https://colab.research.google.com/) : this is
    Google's version of an online notebook, which runs directly in
    Google's cloud.  This starts up almost instantly.

