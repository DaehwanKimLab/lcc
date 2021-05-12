

def Write_Class_Constant_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FConstant():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Variable_("self.NA", 6.022141527E23) # Avogadro's Number

            # Define accessory variables - TF version
            Writer.InitArrayWithOne('OneInt', 1, 'int32')
            Writer.InitArrayWithOne('OneFlt', 1)
            Writer.InitArrayWithZero('ZeroInt', 1, 'int32')
            Writer.InitArrayWithZero('ZeroFlt', 1)
            Writer.BlankLine()

