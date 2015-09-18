'''
Created on 17 Sep 2015

@author: alex
'''
from phe.mapping import Mapper

class BWAMapper(Mapper):
    '''
    classdocs
    '''

    _default_options = ""


    def __init__(self, cmd_options=None):
        '''
        Constructor
        '''

        # Call the
        if not cmd_options:
            cmd_options = self._default_options

        super(BWAMapper, self).__init__(cmd_options=cmd_options)

    def make_sam(self):
        return None


    def get_info(self):
        return None
