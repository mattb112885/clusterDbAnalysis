#!/usr/bin/env python

'''

This file contains basic classes for GUIs based on easygui.
The classes in this file do not depend on any other ITEP scripts 
so if you want to use it for other functions you can as well.

'''

import easygui

import operator
import os
import tempfile

''' Class for exception raised if user cancels an operation.'''
class UserCancelError(Exception):
    pass

''' 
This is the base class for GUI errors displaying error messages.
You can inherit from this class if you want to give a more descriptive name to your errors.
'''
class GuiError(Exception):
    def __init__(self, errormsg):
        msg = "The program encountered the following error:\n\n%s\n\nPress OK to terminate the program.\n" %(errormsg)
        easygui.msgbox(msg=msg)

'''
Error thrown if a gene specified by a user does not exist.
'''
class NoGeneError(GuiError):
    pass

''' This class contains functions that are common to various UI scripts in ITEP '''
class GuiBase:
    def __init__(self):
        # No-op. This class has no attributes.
        pass
    def _createTemporaryFile(self, delete=True):
        ''' Creates a named temporary file and returns the handle to the file and the file's name. 
        
        If delete is specified as True, creates a file that automatically deletes itself when it passes
        out of scope. If delete is specified as False the user will need to delete it manually to make
        it go away (but doesn't need to worry about scope).'''
        f = tempfile.NamedTemporaryFile(delete=delete)
        fname = f.name
        return (f, fname)
    def _save_text(self, text, filename):
        ''' Write text to file filename without manipulating. 
        It is assumed that the calling function already checked for existence. '''
        fid = open(filename, "w")
        fid.write(text)
        fid.close()
        return
    def _success_dialog(self, filename):
        ''' Dialog box for successful saving of a file.'''
        easygui.msgbox(msg = "Successfully saved results to file %s" %(filename) )
        return
    def _get_directory(self):
        ''' Dialog box for getting a directory to which to save files. '''
        return easygui.diropenbox(msg="Where would you like to save the results files?", title="Choose a directory", default=None)
    def _get_file_name(self, extension="txt"):
        '''
        
        This returns the name of a file the user wants to save as
        or None if the user cancels or does not name a file.

        '''
        filename = easygui.filesavebox(msg = "Where do you want to save the file (extension %s will automatically be added)?" %(extension))
        # User cancelled.
        if filename is None:
            return None
        filename = filename + "." + extension
        # Check for file existence and ask if it is OK to overwrite the file.
        if os.path.exists(filename):
            ok_to_overwrite = easygui.buttonbox(msg="File %s already exists. Overwrite?" %(filename), choices = ("No", "Yes") )
            if ok_to_overwrite == "Yes":
                return filename
            else:
                return None
        else:
            return filename
    def _save_file_dialogs(self, extension = "txt"):
        ''' Dialogs asking users to save file, sanity checks for existence of file, etc.

        extension is the file's extension (txt, png, etc). It will automatically be added
        to the filename specified by the user.

        The function returns the filename specified by the user if one is specified or
        None if the user cancels the operation for any reason.
        '''
        # If user cancels it defaults to the FIRST choice. We want default to be NO so I reverse the default of choices here. 
        saveornot = easygui.buttonbox(msg="Do you want to save results to a file?", choices = ("No", "Yes") )
        if saveornot == "Yes":
            filename = self._get_file_name(extension=extension)
            return filename
        else:
            return None
    def _print_readable_table(self, rows, header=True, separator = '|'):
        '''Print a readable table from an array of arrays ("rows").

        Returns a string that can be printed to give you a pretty table.

        Specify header=True to print a header separated from the rest of the table by "-------".
        Separator is the character used to delineate columns. '''
        finaltext = ''
        # What is the maximum length of each column?
        numcols = len(rows[0])
        lens = []
        for i in range(numcols):
            col = map(operator.itemgetter(i), rows)
            maxlen = max(map(len, col))
            lens.append(maxlen)
        # Print the header row if there is one.
        if header:
            headers = rows[0]
            header_elements = []
            sep_elements = []
            for i in range(numcols):
                diff = lens[i] - len(headers[i])
                header_elements.append( ' ' + headers[i] + ' '*diff + ' ')
                sep_elements.append( '-'*( lens[i] + 2 ) )
            finaltext += separator.join(header_elements) + "\n" + separator.join(sep_elements) + "\n"
            rows = rows[1:]
        # Print the rest of the rows.
        formats = []
        for line in rows:
            row_elements = []
            for i in range(numcols):
                diff = lens[i] - len(line[i])
                row_elements.append( ' ' + line[i] + ' '*diff + ' ')
            finaltext += separator.join(row_elements) + '\n'
        return finaltext
