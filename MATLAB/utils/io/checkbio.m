function checkbio()
bio = license('test', 'bioinformatics_toolbox');
fprintf("Test for Bioinformatics Toolbox: %d\n", bio);
assert(logical(bio),"Please Install Bioinformatics Toolbox: https://www.mathworks.com/products/bioinfo.html")
end

