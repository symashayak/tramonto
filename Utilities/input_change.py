
# code to modify dft_input.dat files

# assumes the new style of input file is in dft_input_template
# assumes this file is in the /Examples directory


import string
import os

dirs = os.listdir('.')      # find subdirectories
for dir in dirs :
    if os.path.isdir(dir) : # if this is a directory, look inside
        print 'fixing directory %s\n' % dir,
        nofile = 0

        # if the file dir/dft_input.dat exists, then modify it
        oldfilename = dir + '/dft_input.dat'
        try :
            fileold = open(oldfilename,'r')
        except IOError :
            nofile = 1

        if nofile == 0 :
        
            filenew = open('dft_input_template', 'r')
            modfilename = dir + '/input_mod'
            file_mod = open(modfilename,'w')

            while 1:            # loop over lines in both files; they can have different numbers of lines
                while 1 :
                    lineold = fileold.readline()
                    if not lineold : break
                    elif lineold[0] == '@' :
                        break
                while 1:
                    linenew = filenew.readline()
                    if not linenew: break
                    elif linenew[0] == '@' :
                        break
                    else :
                        file_mod.write(linenew)         # echo the new line if it's not an "@" line
             
                    
                if not linenew : break
                if not lineold : break

                count = 0
                flag = 0
                if lineold[0] == '@':                   # process each line; keep the numbers from fileold and insert the text from filenew
                        cols_old = lineold.split()      # split into space-delimited "words"
                        cols_new = linenew.split()
                     #   print cols_old
                        file_mod.write('@ ')
                        for item in cols_old :
                            if item == '@' :
                                continue
                            if not (item[0] in string.letters) :    # if it's not a letter, it must be a number; keep it
                      #          print item
                                item = item + ' '
                                file_mod.write(item)
                            elif item == 'n/a' :                     # special case
                                item = item + ' '
                                file_mod.write(item)
                            else :                                  # once we hit the first non-number, switch to the filenew
                                flag = 1

                            if flag :
                                break

                        if flag :
                            for item in cols_new :
                                if item == '@' :
                                    count += 1
                                elif not (item[0] in string.letters) :
                                    count += 1
                                else :
                                    break
                            newtext = cols_new[count:]
                            file_mod.write('\t')
                            for item in newtext :
                                item = item + ' '
                                file_mod.write(item)

                        file_mod.write('\n')

            fileold.close()
            filenew.close()
            file_mod.close()

            # change file names; keep old version just in case
            mvold = dir + '/input_old.dat'
            os.rename(oldfilename,mvold)
            os.rename(modfilename,oldfilename)

print '\ndone\n'


