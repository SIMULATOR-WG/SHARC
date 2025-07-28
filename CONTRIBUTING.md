# Contributing

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

## Types of Contributions

### Report Bugs

Report bugs at https://github.com/Radio-Spectrum/SHARC/issues

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

### Fix Bugs
Look through the GitHub issues for bugs. Anything tagged with "bug"
and "help wanted" is open to whoever wants to implement it.

### Implement Features
Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

### Write Documentation
SHARC could always use more documentation, whether as part of the
official SHARC docs, in docstrings, or even on the web in blog posts,
articles, and such.

When implementing a new feature such as a new propagation model or a new antenna model,
add the documentation to the Wiki.

### Submit Feedback
The best way to send feedback is to file an issue at https://github.com/Radio-Spectrum/SHARC/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

## Branching model

The branching model is based on this excelent post written by Vincent Driessen.

A Successful Git Branching Model - https://nvie.com/posts/a-successful-git-branching-model/

### The main branches
The central repo holds two main branches with an infinite lifetime:

* `master` - main branch where the source code of HEAD always reflects a production-ready state.
* `development` - the main branch where the source code of HEAD always reflects a state with the latest 
delivered development changes for the next release

When the source code in the `development` branch reaches a stable point and is ready to be released, 
all of the changes should be merged back into `master` and then tagged with a release number.

### supporting branches
The supporting porting branches to aid parallel development between team members, ease tracking of features, 
prepare for production releases and to assist in quickly fixing live production problems. 
Unlike the main branches, these branches always have a limited life time, since they will be removed eventually.

* `feature` branches
* `release` branches
* `hotfix` branches

### feature branches
Feature branches are used to develop new features. For example, an addidion of a new antenna model or
a new propagation model. 

May branch off from:
    `develop`
Must merge back into:
    `develop`
Branch naming convention:
    anything except `master`, `development`, `release-*`, or `hotfix-*`. For example, `feature-propagation-619`.

### release branches
Release branches support preparation of a new production release. They allow for minor bug fixes and 
preparing meta-data for a release (version number, build dates, etc.)

May branch off from:
    `develop`
Must merge back into:
    `develop` and `master``
Branch naming convention:
    `release-*`

### hotfix branches
Used when a critical bug in a production version must be resolved immediately.

May branch off from:
    `master`
Must merge back into:
    `develop` and `master`
Branch naming convention:
    `hotfix-*`

## Languages

### Python

The [PEP8 style guide](https://www.python.org/dev/peps/pep-0008/) is authoritative.

#### Type annotations
All new code should be fully type-annotated.
For reference, please look at this [type hints cheat sheet for Python 3](https://mypy.readthedocs.io/en/stable/cheat_sheet_py3.html).

#### Documentation
* Document all public functions and keep those docs up to date when you make changes
* We use [Google style docstrings](https://google.github.io/styleguide/pyguide.html#383-functions-and-methods) in our codebase (if you find a docstring missing or in the wrong format you're welcome to fix it.)
    * For VSCode users, [Python Docstring Generator](https://marketplace.visualstudio.com/items?itemName=njpwerner.autodocstring) plugin is recommended
Example:

```
def foo(arg1: str) -> int:
    """Returns the length of arg1.

    Args:
        arg1 (str): string to calculate the length of

    Returns: the length of the provided parameter
    """
    return len(arg1)
```


## Get Started!
Ready to contribute? Here's how to set up `sharc` for local development.

1. Fork the `sharc` repo on GitHub.
2. Clone your fork locally:

    `$ git clone git@github.com:your_name_here/sharc.git`

    `$ cd SHARC/`

3. Install any version of python 3 and the virutalenv module in your system (SHARC has been tested from versions 3.8 and above).
4. Install your local copy into a virtualenv.
    
    `$ cd sharc`

    `$ python3 -m venv .venv`

    `$ source .venv/bin/activate`

    You shall see (.venv) in the begining of your command prompt indicating that the virutalenv has been activated.
    Now, instalal the dependencies for development.
    `$ pip install -r requirements.txt`

    Install sharc on your local enviroment.
    Run from the source code directory root:
    `$ pip install -e .`

4. Create a branch for local development::

    `$ git checkout -b name-of-your-bugfix-or-feature`

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the tests, including testing other Python versions with tox:

    `$ flake8 sharc tests`

    `$ python setup.py test or py.test`


   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub::

    `$ git add .`

    `$ git commit -m "Your detailed description of your changes."`
    
    `$ git push origin name-of-your-bugfix-or-feature`

7. Submit a pull request through the GitHub website.

# Pull Request Guidelines

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 2.6, 2.7, 3.3, 3.4 and 3.5, and for PyPy. Check
   https://travis-ci.org/edgar-souza/sharc/pull_requests
   and make sure that the tests pass for all supported Python versions.

Tips
----

To run a subset of tests::

$ py.test tests.test_sharc

