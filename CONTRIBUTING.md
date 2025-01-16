# Contributors Guide

**Welcome to the Agate repository!**
We're excited you're here and want to contribute.

We hope that these guidelines make it as easy as possible to get involved.
If you have any questions that aren't discussed below, please let us know by:

* [Opening a GitHub issue](https://github.com/agate-model/Agate.jl/issues/new)

* [Creating a GitHub discussion](https://github.com/agate-model/Agate.jl/discussions/new/choose)

We welcome all contributions from documentation to testing to writing code.
Don't let trying to be perfect get in the way of being good - exciting ideas are more important than perfect pull requests.

## Where to start: issues

The simplest way to contribute to Agate is to create or comment on issues and discussions. Before you open a new issue, please check if any of our [open issues](https://github.com/agate-model/Agate.jl/issues) covers your idea already.

Particularly useful issues for us are bug reports:

* Provide an explicit code snippet --- not just a link --- that reproduces the bug in the latest tagged version of Agate. This is sometimes called the ["minimal working example"](https://en.wikipedia.org/wiki/Minimal_working_example). Reducing bug-producing code to a minimal example can dramatically decrease the time it takes to resolve an issue.

* Paste the _entire_ error received when running the code snippet, even if it's unbelievably long.

* Use triple backticks (e.g., ````` ```some_code; and_some_more_code;``` `````) to enclose code snippets, and other [markdown formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) to make your issue easy and quick to read.

* Report the Agate version, Julia version, machine and any other possibly useful details of the computational environment in which the bug was created.

Discussions are recommended for asking questions about (for example) the user interface, implementation details, science, and life in general.

## Making a change with a pull request

We appreciate all contributions to Agate.
**THANK YOU** for helping us.

We follow the [ColPrac guide](https://github.com/SciML/ColPrac) for collaborative practices.
We ask that new contributors read that guide before submitting a pull request.

The following steps are a guide to help you contribute in a way that will be easy for everyone to review and accept with ease.

### 1. Comment on an [existing issue](https://github.com/agate-model/Agate.jl/issues) or open a new issue referencing your addition

This allows other members of the Agate team to confirm that you aren't overlapping with work that's currently underway and that everyone is on the same page with the goal of the work you're going to carry out.

[This blog](https://www.igvita.com/2011/12/19/dont-push-your-pull-requests/) is a nice explanation of why putting this work in up front is so useful to everyone involved.

### 2. [Fork][github-fork] the [Agate repository][Agate-repo] to your profile

This is now your own unique copy of Agate.
Changes here won't affect anyone else's work, so it's a safe space to explore edits to the code!

Make sure to [keep your fork up to date][github-syncfork] with the master repository, otherwise you can end up with lots of dreaded [merge conflicts][github-mergeconflicts].

### 3. Make the changes you've discussed

Try to keep the changes focused.
If you submit a large amount of work all in one go it will be much more work for whomever is reviewing your pull request.

While making your changes, commit often and write good, detailed commit messages.
[This blog](https://chris.beams.io/posts/git-commit/) explains how to write a good Git commit message and why it matters.
It is also perfectly fine to have a lot of commits - including ones that break code.
A good rule of thumb is to push up to GitHub when you _do_ have passing tests then the continuous integration (CI) has a good chance of passing everything.

If you feel tempted to "branch out" then please make a [new branch][github-branches] and a [new issue][Agate-issues] to go with it. [This blog](https://nvie.com/posts/a-successful-git-branching-model/) details the different Git branching models.

Please do not re-write history!
That is, please do not use the [rebase](https://help.github.com/en/articles/about-git-rebase) command to edit previous commit messages, combine multiple commits into one, or delete or revert commits that are no longer necessary.

### 4. Submit a [pull request][github-pullrequest]

We encourage you to open a pull request as early in your contributing process as possible.
This allows everyone to see what is currently being worked on.
It also provides you, the contributor, feedback in real time from both the community and the continuous integration as you make commits (which will help prevent stuff from breaking).

When you are ready to submit a pull request, make sure the contents of the pull request body do the following:
- Describe the problem you're trying to fix in the pull request, reference any related issues and use keywords fixes/close to automatically close them, if pertinent.
- List changes proposed in the pull request.
- Describe what the reviewer should concentrate their feedback on.

If you have opened the pull request early and know that its contents are not ready for review or to be merged, add "[WIP]" at the start of the pull request title, which stands for "Work in Progress".
When you are happy with it and are happy for it to be merged into the main repository, change the "[WIP]" in the title of the pull request to "[Ready for review]".

A member of the Agate team will then review your changes to confirm that they can be merged into the main repository.
A [review][github-review] will probably consist of a few questions to help clarify the work you've done.
Keep an eye on your GitHub notifications and be prepared to join in that conversation.

You can update your [fork][github-fork] of Agate [repository][Agate-repo] and the pull request will automatically update with those changes.
You don't need to submit a new pull request when you make a change in response to a review.

You can also submit pull requests to other contributors' branches!
Do you see an [open pull request](https://github.com/agate-model/Agate.jl/pulls) that you find interesting and want to contribute to?
Simply make your edits on their files and open a pull request to their branch!

What happens if the continuous integration (CI) fails (for example, if the pull request notifies you that "Some checks were not successful")?
The CI could fail for a number of reasons.
At the bottom of the pull request, where it says whether your build passed or failed, you can click “Details” next to the test, which takes you to the GitHub Actions page.

GitHub has a [nice introduction][github-flow] to the pull request workflow, but please get in touch if you have any questions.

## Setting up your development environment

* Install [Julia](https://julialang.org/) on your system.

* Install `git` on your system if it is not already there (install XCode command line tools on
  a Mac or `git bash` on Windows).

* Login to your GitHub account and make a fork of the
  [Agate repository](https://github.com/agate-model/Agate.jl) by
  clicking the "Fork" button.

* Clone your fork of the Agate repository (in terminal on Mac/Linux or git shell/
  GUI on Windows) in the location you'd like to keep it.
  ```
  git clone https://github.com/your-user-name/Agate.jl.git
  ```

* Create the development environment by opening Julia via `julia --project` then
  typing in `] instantiate`. This will install all the dependencies in the Project.toml
  file.

* You can test to make sure Agate works by typing in `] test`. Doing so will run all
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

_These Contributing Guidelines have been adapted directly from [The Turing Way](https://github.com/the-turing-way/the-turing-way/blob/main/CONTRIBUTING.md) and [OceanBioME](https://github.com/OceanBioME/OceanBioME.jl/blob/main/docs/src/contributing.md)_


[Agate-repo]: https://github.com/agate-model/Agate.jl
[Agate-issues]: https://github.com/agate-model/Agate.jl/issues
[git]: https://git-scm.com
[github]: https://github.com
[github-branches]: https://help.github.com/articles/creating-and-deleting-branches-within-your-repository
[github-fork]: https://help.github.com/articles/fork-a-repo
[github-flow]: https://guides.github.com/introduction/flow
[github-mergeconflicts]: https://help.github.com/articles/about-merge-conflicts
[github-pullrequest]: https://help.github.com/articles/creating-a-pull-request
[github-review]: https://help.github.com/articles/about-pull-request-reviews
[github-syncfork]: https://help.github.com/articles/syncing-a-fork