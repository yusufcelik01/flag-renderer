# flag-renderer
A simple waving flag rendering openGL program written in cpp for ceng469 computer graphics 2 course. 

# usage
```
make flag
./flag [image_file]
```
should create an window containing an animated flag continuously. If no input image file is given the program searches for metu_flag.png in the current directory and uses it for texture mapping.


# controls
pressing **key W** and **key S** increases or decrases sampling rate respectively. Low sample rates would show sharp edges and high sample rates may slow down the program.

Pressing **key E** increases the number of bezier patches in the flag and **key D** decreases the number of beziers patches.

# sample output
![snapshot of programs output](https://github.com/yusufcelik01/flag-renderer/blob/main/render_output.png)
