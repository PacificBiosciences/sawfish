# Sawfish methods documentation

The sawfish methods documentation is markdown with inline latex equations. It is designed for rendering through pandoc.
With pandoc and latex installed, a pdf methods document can be generated as follows:

    pandoc methods.md -f markdown -t pdf -o methods.pdf

If these tools are not installed, one approach is to run the above command through a docker pandoc image as follows:

    docker run --rm -v "$(pwd)":/data -w /data -u $(id -u):$(id -g) pandoc/latex methods.md -f markdown -t pdf -o methods.pdf
