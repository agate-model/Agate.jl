# Contributors Guide

**Welcome to the Agate repository!**
We're excited you're here and want to contribute.

We hope that these guidelines make it as easy as possible to get involved.
If you have any questions that aren't discussed below, please let us know by:

  - [Opening a GitHub issue](https://github.com/agate-model/Agate.jl/issues/new)

  - [Creating a GitHub discussion](https://github.com/agate-model/Agate.jl/discussions/new/choose)

We welcome all contributions from documentation to testing to writing code.
Don't let trying to be perfect get in the way of being good - exciting ideas are more important than perfect pull requests.

## Where to start: issues

The simplest way to contribute to Agate is to create or comment on issues and discussions. Before you open a new issue, please check if any of our [open issues](https://github.com/agate-model/Agate.jl/issues) covers your idea already.

Particularly useful issues for us are bug reports:

  - Provide an explicit code snippet --- not just a link --- that reproduces the bug in the latest tagged version of Agate. This is sometimes called the ["minimal working example"](https://en.wikipedia.org/wiki/Minimal_working_example). Reducing bug-producing code to a minimal example can dramatically decrease the time it takes to resolve an issue.

  - Paste the _entire_ error received when running the code snippet, even if it's unbelievably long.
  - Use triple backticks (e.g., ````some_code; and_some_more_code;````) to enclose code snippets, and other [markdown formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) to make your issue easy and quick to read.
  - Report the Agate version, Julia version, machine and any other possibly useful details of the computational environment in which the bug was created.

Discussions are recommended for asking questions about (for example) the user interface, implementation details, science, and life in general.

## Making a change with a pull request

We appreciate all contributions to Agate.
**THANK YOU** for helping us.

We follow the contributing guidelines developed by [The Turing Way](https://github.com/the-turing-way/the-turing-way/blob/main/CONTRIBUTING.md#making-a-change-with-a-pull-request). Please read the full guide, especially if you are new to contributing to open source projects. In short, there are four steps to adding changes to this repository:

 1. **Capture the change in an issue**: Comment on an [existing issue](https://github.com/agate-model/Agate.jl/issues) or open a new issue referencing your addition.
 2. **Fork the Repository**: [Fork the Agate repository](https://github.com/agate-model/Agate.jl/fork).
 3. **Make Changes**: Commit often and write good, detailed commit messages (see [this blog](https://chris.beams.io/posts/git-commit/)).
 4. **Open a Pull Request**: Ensure you describe the changes made and what the reviewer should focus on as well as any additional details.

## Setting up your development environment

  - Install [Julia](https://julialang.org/) on your system.

  - Install `git` on your system if it is not already there (install XCode command line tools on
    a Mac or `git bash` on Windows).
  - Login to your GitHub account and make a fork of the
    [Agate repository](https://github.com/agate-model/Agate.jl) by
    clicking the "Fork" button.
  - Clone your fork of the Agate repository (in terminal on Mac/Linux or git shell/
    GUI on Windows) in the location you'd like to keep it.
    
    ```
    git clone https://github.com/<your-user-name>/Agate.jl.git
    ```
  - Create the development environment by opening Julia via `julia --project` then
    typing in `] instantiate`. This will install all the dependencies in the Project.toml
    file.
  - You can test to make sure Agate works by typing in `] test`. Doing so will run all
    the tests (and this can take a while).

We follow the [Blue](https://github.com/JuliaDiff/BlueStyle) style guide for Julia. To automatically format all Julia files in the project you can use the JuliaFormatter. Once you have installed it (`add JuliaFormatter`) run:

```Julia
using JuliaFormatter

format(".")
```

Your development environment is now ready!

## Documentation

Now that you've made your awesome contribution, it's time to tell the world how to use it.
Writing documentation strings is really important to make sure others use your functionality
properly. Didn't write new functions? That's fine, but be sure that the documentation for
the code you touched is still in great shape. It is not uncommon to find some strange wording
or clarification that you can take care of while you are here.

You can preview how the Documentation will look like after merging by building the documentation
locally. From the main directory of your local repository call

```
JULIA_DEBUG=Documenter julia --project=docs/ docs/make.jl
```

and then open `docs/build/index.html` in your favorite browser. Providing the environment variable
`JULIA_DEBUG=Documenter` will provide with more information in the documentation build process and
thus help figuring out a potential bug.

## Credits

_These Contributing Guidelines have been adapted directly from [The Turing Way](https://github.com/the-turing-way/the-turing-way/blob/main/CONTRIBUTING.md) and [OceanBioME](https://github.com/OceanBioME/OceanBioME.jl/blob/main/docs/src/contributing.md)._
