# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 18:57:33 2012

@author: mag
"""

class adder:

    result = 0

    def __init__( self, number1, number2 ):
        self.result = int( number1 ) + int( number2 )

    def giveResult( self ):
        return str(self.result)

    endIt = False
    while ( endIt == False ):
        print "Please input two intergers you wish to add: "
        number1 = raw_input( "Enter the first number: " )
        number2 = raw_input( "Enter the second number: " )
        try:
            thistime = adder( number1, number2 )
        except ValueError:
            print "Sorry, one of your values was not a valid integer."
        continue
        print "Your result is: " + thistime.giveResult()
        goagain = raw_input( "Do you want to eXit or go again? ('X' to eXit, anything else to continue): " )
        if ( goagain == "x" or goagain == "X" ):
            endIt = True