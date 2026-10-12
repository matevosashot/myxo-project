
# Matheamtica
do not edit any mathematica file (nb, wls, m), unless explicitly asked to do so.

Only prompt the changes and replacements.

I allow you to create and run wolfram scripts in ./claude_experiments. If you are unsure about some syntax in mathematica, you may write a script and run it using wolframscript.

module load mathematica/14.10 to enable mathematica.

In order to export an image with wolframscript, unset DISPLAY environment variable.

# Latex

## Highlight the changes

Use \upd{...} for smaller changes or \begin{updenv}... \end{updenv} for larger blocks of text to highlight the changes.


## Prompts in latex
You do this, when I specifically ask to address the prompts. If you are unsure, ask.

In .tex files, you will find commands in the form (can be commented out)
\prompt{....}

You should read this prompt and address it, even if it is commented out.

Do not address \todo commands, unless I specifically ask to; leave them as they are.

Address each comment thoroughly. Each prompt can affect multiple parts of the manuscript, not only local changes.


So for each prompt:
* understand what is it about.
* Ask questions if it is not clear how to address it, with suggestions
* plan the change
* apply the change
* remove that prompt from the tex
* update the memory if needed

Highlight the changes (e.g. text color dark green) so I can see what you changed.

If there are existing highlights, ask (1) remove the previous highlights, (2) use different colors (2) leave them as they are.