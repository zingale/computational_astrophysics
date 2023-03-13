# github

Github provides a web-based way to interact with git repositories.  At
its heart it hosts a bare repo that we can push-pull to/from, but it
also provides features to make it easier for users to request their
changes be incorporated.


## Creating a repository on github

Let's start by creating a new git repository using github's web interface.  Start
on your github home page and click on the "+" icon and select "New repository":

```{image} github-new.png
:align: center
```

Now we give the repository a name.  Let's use our initials, followed
by `_class_repo`, so for me, it will be `mz_class_repo`.

Make sure that it defaults the repo to be public, which means anyone on the internet
can see the contents&mdash;that's what we want.

Finally, check the box to add a `README` file&mdash;this means that our repository will
not be empty initially.

```{image} github-create.png
:align: center
```

Our project is now found at: *https://github.com/username/reponame/*,
where *username* is your Github username and *reponame* is the name of
the repository you just created.


## SSH interlude

Github works best is we communicate via [secure
shell](https://en.wikipedia.org/wiki/Secure_Shell) or *SSH*.

There is some nice documentation describing key pairs in the Software
Carpentry [Create an SSH key pair](https://swcarpentry.github.io/git-novice/07-github/index.html#3-ssh-background-and-setup)
section.

Here's how we will set things up:

1. A the bash prompt generate a new key pair:

   ```bash
   ssh-keygen -t ed25519
   ```

   The `-t` option picks a secure encryption method.

   It will ask you for a passpharse&mdash;just hit "Enter" to keep it
   empty (if other people had access to your account, the you would
   want to pick a passphrase).

   If you do

   ```bash
   ls -l ~/.ssh
   ```

   you'll see 2 files: `id_ed25519` and `id_ed25519.pub` this is
   the private and public key for encryption.

   ```{caution}
   Never share your private key (`id_ed25519`) with anyone.

   Our public key (`id_ed25519.pub`) is meant to be public, and
   we can give it to places we want to communicate with, like github
   ```

2. Go to you Github profile SSH keys settings: https://github.com/settings/keys

   Click on the *New SSH key* button and:

   * give a *title* which is descriptive of the machine you are using, like
     ``MathLab``

   * copy and paste the contents of `id_ed25519.pub` into the *key*
     text box.  You can see the contents by doing:

     ```bash
     cat ~/.ssh/id_ed25519.pub
     ```
     
   * Click on ``Add SSH key``

3. Test things out:

   ```bash
   ssh -T git@github.com
   ```

   It will ask you if we want to save the *fingerprint*&mdash;say "yes", and then
   it should report:

   ```
   Hi zingale! You've successfully authenticated, but GitHub does
   not provide shell access.
   ```

That means everything is working.


## Working remotely

Now we can git clone this repo.  From the github project page, click on the
*code* button.

```{image} github-clone.png
:align: center
```

Copy the string in the text box there and then on your command line clone
the repo as:

```bash
git clone git@github.com:zingale/mz_class_repo.git
```

(replacing my repo and username with your own).

Now we can go into our repo and look around.  Notice that there is a
`.git/` directory.  Also look at the remotes:

```bash
git remote -v
```

```
origin	git@github.com:zingale/mz_class_repo.git (fetch)
origin	git@github.com:zingale/mz_class_repo.git (push)
```

This is just like the example or remotes we did previously, except now
github is acting as our remote.

This means that we call push to github and pull from there.

As a single user, this will allow you to develop from any computer
and keep the code base in sync across all of them.

If the project has multiple developers, this can be where all of the
developers sync up their projects.


### `README.md` is special

The web interface that github provides to our repo has a number of features.

First, the `README.md` file is always displayed on the main project
page.  This is where you can put descriptions of what your project is,
how people can contribute, even share the status of testing and
documentation builds (we'll talk about those later in class).

This file is in github-flavored [Markdown
format](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax)
(that's what the `.md` extension signifies).  Markdown allows you to
do basic formatting.

Here's an example of what you can do in a `README.md` from one of my
projects: https://github.com/pynucastro/pynucastro


## Pull requests

Github allows you to give permissions to users to contribute to a
repository (read, write, or admin access).  But the best type of workflow
is one where users don't push directly to the git repo.  Instead it is based
around pulls.

How do we contribute to a project that we don't own?

Here's a github *organization* for our class: https://github.com/phy504-sbu

and here's a simple repo in this organization: https://github.com/phy504-sbu/test-repo

An organization is meant to be used by a collection of developers who
can all have different access permissions.  It provides tools for
managing who can do different things to the repos under its control.

Let's clone this repo:

.. prompt:: bash

   git clone git@github.com:phy504-sbu/test-repo.git
   cd test-repo

Now, let's each try to add a file of the form *username.txt* containing
your full name.  Ex:

.. prompt:: bash

   echo Michael Zingale > zingale.txt
   git add zingale.txt
   git commit

Now try to push it to the repo we clone:

.. prompt:: bash

   git push

what happened?

The issue is that you don't have *write* permission to that repo,
since I own it.  So you are denied access.

This is okay.  The workflow that github emphasizes is one based around
*pulls* not *pushes*, so let's see how we do that.

First, we need to *fork* the repo -- this creates a clone under our
control that we can do with as we please.  Click on the "fork" button.

.. image:: github-fork.png
   :align: center

It may ask you where you want the fork to life -- you want it to live
under your profile / username.

This will bring you to a new github page, displaying the fork, with a
URL that should look something like: https://github.com/zingale/test-repo

Now click on the *code* button and copy the SSH location.

We want to add this fork as a new remote:

.. prompt:: bash

   git remote add myfork git@github.com:zingale/test-repo.git

(again, make sure you replace that with the link to your repo).

Now you can do:

.. prompt:: bash

   git push myfork

If you reload your github page, you should see your change there.

Now we can do all *pull-request*.  Select "pull requests"

.. image:: github-pr.png
   :align: center

Then click on the "New pull request" button, and you'll see something like:

.. image:: github-pr2.png

This is showing that you are asking to merge the changes in your fork into the
class ``test-repo`` repository.

Click on *create pull request*, type in a comment about what this does, and then click
on the *create pull request* button again.

Now it is out of your hands.

The owner of the class repo (me) will get a notification that you want
to incorporate your changes into the class repo, and I can merge them
via the github web tools.


The overall workflow that we did: fork, push to our fork, issue a PR, looks like:

.. image:: github-workflow.png
   :align: center
   :width: 80%




Our class notes github
======================

Let's take a tour of our class notes on github: https://github.com/zingale/phy504

There are a lot of other features that github provides that we will explore later in the semester, including:

* github actions : automating some workflows (like testing) on our code

* github pages : building and hosting web pages for our project


