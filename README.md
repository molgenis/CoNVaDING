 ![image](https://travis-ci.org/duartemolha/CoNVaDING_reload.svg?branch=master)

# CoNVaDING_reload

This is a fork of the original project CoNVaDING
This project is the official one and I make no claims to be a replacement of it.

If you want to follow the official project please go to: 

https://molgenis.gitbook.io/convading/

# History

This fork was created by Duarte Molha to implement features not provided by the original project.

You can read more about how it differs from the original software here:

https://github.com/duartemolha/CoNVaDING_reload/wiki/CoNVaDING_reload

Sequentially again forked by mmterpstra who added an amplicon mode to CountFromBam, implemented the perl build system and misc bugfixes.

# Install

How to install on ubuntu 18.04 to your home:

```
# install bedtools / samtools 
sudo apt install bedtools samtools

# install convading
cpanm --installdeps --notest .
perl Makefile.PL PREFIX=~/
make
make install
```


# Developing

```
# change depending on your preferences

# preferably fork first
git clone https://github.com/duartemolha/CoNVaDING_reload.git
git branch feature-tag
git checkout feature-tag
# edit some files
# test by using the following commands

(export PERL5LIB="blib/lib/:$PERL5LIB";perl Makefile.PL && make && make test)

# add, commit and push changed files
# create pull request on github
```

