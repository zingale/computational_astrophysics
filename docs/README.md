# Building the course notes

To build, as-is, without running the notebooks, do:

```
make SKIP_EXECUTE=TRUE html
```

if you want to execute the notebooks (this can take ~10 minutes), then
do:

```
make html
```

The site will then be in `build/html/`.  You can open it with

```
firefox build/html/index.html
```

or, if you want to mimic a real webserver, as:

```
python3 -m http.server 8000 --bind 127.0.0.1 --directory build/html
```

and then point your browser to http://localhost:8000
