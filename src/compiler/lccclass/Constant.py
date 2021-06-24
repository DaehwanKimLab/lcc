
# Comp is a short hand for CompilerData
def Write_Constant(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FConstant():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Variable_("self.NA", 6.022141527E23) # Avogadro's Number

            # Define accessory variables - TF version
            Writer.InitOnes_('OneInt', 1, 'int32')
            Writer.InitOnes_('OneFlt', 1)
            Writer.InitZeros('ZeroInt', 1, 'int32')
            Writer.InitZeros('ZeroFlt', 1)
            Writer.BlankLine()

