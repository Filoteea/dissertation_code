# Filo's Coding Etiquette :nerd_face:

## Script structure  :page_with_curl:
### The script will always be divided in sections (*using ---- at the end of the comment line*)

  - **Introduction**: Author statement (what does this script do?), author(s) names, contact details and date. 

  - **Libraries**: Keep all the libraries within this section. Add comments to each library explaining what are you using each package for. 
 
  - **Functions**: Define all the functions used throughout the code, add comments explaining what they do and mentioning who wrote them if it's not you. When you open RStudio, remember to run this section.

  - **Set the working directory**: Will show the filepath to the WD. The WD will contain all the files needed for the project structured within the folders: scripts, inputs, outputs, figures, raw_data, research_references.

  - **Import the data**: All the data used for a project will be stored in a folder in the working directory. Add comments describing the data you are using (format, variable type, size) and explain what you will use it for.

  - **The different sections of your analysis**: Create a section for each part of the analysis in a cohesive logical order. 

  - **The outputs of your analysis**: Save figures or other outputs in the designated folder. Any images will be saved in both .pdf and .png formats.  

## Giving names :baby:

 - **Variables**: Only use nouns separated with underscores, written in lower case (e.g. *methane_lifetime*).
  
 - **Functions**: Use verbs, words are separated with dots, written in lower case (e.g. *calculate.methane.lifetime*). 
  
 - **Files**: The names have to explain the purpose of the file, mention the date if project related to annual/ motnhly etc. data. Separated with underscore and can use both lower case and upper case (e.g. *CH4_lifetime_Oct_2021.R; CH4_emissions_Sept_2021.csv*).

## Code syntax rules :warning:

### Spacing 

 :heavy_check_mark: Place spaces around all infix operators (=, +, -, <-, etc.). The same rule applies when using = in function calls. 
 
 :heavy_check_mark: Always put a space after a comma, and never before.

 :heavy_check_mark: Exceptions:  ":" and "::" don’t need spaces around them; you shouldn't add spaces when defining coordinate systems in spatial objects.
             
 :heavy_check_mark: Don’t place a space before left parentheses, except in a function call.
 - *Example*: 

``` {r}
 if (debug) do(x)      
 plot(x, y)
```     

:heavy_check_mark: Always place four spaces between code and comment.

### Curly braces 

 :heavy_check_mark: An opening curly brace should never go on its own line and should always be followed by a new line.
 
 :heavy_check_mark: A closing curly brace should always go on its own line, unless it’s followed by else. Always indent the code inside curly braces.
 - *Example*:   
```{r}
if (y < 0 && debug) {      
  message("Y is negative")       
}    
    
if (y == 0) {      
  log(x)     
} else {     
  y ^ x    
}
```      

### Other 
:heavy_check_mark: **Indentation**: If a command runs over multiple lines, indent the second line to where the definition starts.     
:heavy_check_mark: Always use " instead of '.




