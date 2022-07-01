## Download the fastq files for workshop

Goto https://cloudstor.aarnet.edu.au/plus/s/ZAwa5VsRZYjqpai
- select all the files 
- click on "Download"

The files downlaod to your laptop under "subsamples.tar" 

## Let's transfer them to the server provided for workshop
On your laptop, use WinSCP https://winscp.net/eng/index.php to transfer file

If you are using command line, on your laptop 

    scp –r subsample.tar username@115.146.84.253:~/.
  
## Confirm files were transfered to the server 
Log in to the server 
  - use ssh command line 
  
        ssh user@115.146.84.253
    
  - Use putty or Mobaterm to login instead
  
Once logged in, type the command 

    ls 
  
This should list the folder "subsample.tar" 


## Decompress the file
Type in the below command 

    tar -xvzf subsample.tar
    ls
  
The decompressed folder will be labelled "subsample". 


Yay \0/, all done with this tutorial !
