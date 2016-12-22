# LeachingSimulations
Contaminant Transport in Microcosms

This repository collects progress on leaching models written in python.  

To contribute follow these steps: http://rogerdudler.github.io/git-guide/

A. Download the project to your computer

1. Fork the project to your own GitHub account (see fork button on top right corner of GitHub page)
2. In the terminal change to the directory where you want to store the project
3. Clone it by typing in the command window: "git clone" + "http address of the repository" (note: must have Git installed)
4. Make the changes needed in your computer.

B. Commit Changes

1. After you've made changes to the files in local directory (i.e. your computer), check the status on your terminal by typing :"git status"
2. Add the files you edited: "git add <filename>" or to add all files: "git add *"
3. To actually commit these changes use: " git commit -m "Commit message" ". Now the file is committed to the HEAD, but not in your remote repository yet.

C. Push Changes

1. To send those changes to your remote repository, execute: "git push origin master". Change master to whatever branch you want to push your changes to.

D. Branching

1. Make new branch locally: git checkout -b new_branch_name
2. Switch back to master: git checkout master

3. Delete the branch locally: git branch -d new_branch_name
4. Delete the branch remotely: git push origin --delete new_branch_name
or:
5. Make it available remotely: git push origin new_branch_name
