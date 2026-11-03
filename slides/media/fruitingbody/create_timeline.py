import os

files = [file for file in os.listdir(".")
         if os.path.splitext(file)[1]==".png"]

with open("timeline.txt", "w") as f:
    for i,file in enumerate(files):
        f.write(f"::{i}\n")
