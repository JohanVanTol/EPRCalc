object OptionsDialog: TOptionsDialog
  Left = 192
  Top = 107
  Caption = 'OptionsDialog'
  ClientHeight = 343
  ClientWidth = 330
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBox1: TGroupBox
    Left = 8
    Top = 8
    Width = 289
    Height = 113
    Caption = 'Sigma square calculation'
    TabOrder = 1
    object Label1: TLabel
      Left = 16
      Top = 16
      Width = 169
      Height = 13
      Caption = 'Number of Points for  noise estimate'
    end
    object Label9: TLabel
      Left = 16
      Top = 74
      Width = 148
      Height = 13
      Caption = 'Data Column for Noise estimate'
    end
    object NoisePointsEdit: TEdit
      Left = 216
      Top = 13
      Width = 49
      Height = 21
      TabOrder = 0
      Text = '30'
    end
    object NoiseColumnComboBox: TComboBox
      Left = 216
      Top = 66
      Width = 49
      Height = 21
      Style = csDropDownList
      ItemHeight = 13
      TabOrder = 1
    end
    object FromFirstRadioButton: TRadioButton
      Left = 16
      Top = 40
      Width = 113
      Height = 17
      Caption = 'From First'
      TabOrder = 2
    end
    object FromLastRadioButton: TRadioButton
      Left = 144
      Top = 40
      Width = 113
      Height = 17
      Caption = 'From Last'
      TabOrder = 3
    end
  end
  object BitBtn1: TBitBtn
    Left = 239
    Top = 307
    Width = 75
    Height = 25
    TabOrder = 0
    Kind = bkOK
  end
  object GroupBox2: TGroupBox
    Left = 8
    Top = 147
    Width = 225
    Height = 185
    Caption = 'Optimization parameters'
    TabOrder = 2
    object Label2: TLabel
      Left = 32
      Top = 32
      Width = 92
      Height = 13
      Caption = 'Change of g-values'
    end
    object Label3: TLabel
      Left = 16
      Top = 56
      Width = 106
      Height = 13
      Caption = 'Change of ZFS values'
    end
    object Label4: TLabel
      Left = 192
      Top = 56
      Width = 21
      Height = 13
      Caption = 'GHz'
    end
    object Label5: TLabel
      Left = 32
      Top = 80
      Width = 93
      Height = 13
      Caption = 'Change of A-values'
    end
    object Label6: TLabel
      Left = 192
      Top = 80
      Width = 22
      Height = 13
      Caption = 'MHz'
    end
    object Label7: TLabel
      Left = 192
      Top = 104
      Width = 15
      Height = 13
      Caption = 'mT'
    end
    object Label8: TLabel
      Left = 32
      Top = 104
      Width = 77
      Height = 13
      Caption = 'Change of width'
    end
    object gChangeEdit: TEdit
      Left = 136
      Top = 24
      Width = 73
      Height = 21
      TabOrder = 0
      Text = '0.00003'
    end
    object ZFSChangeEdit: TEdit
      Left = 136
      Top = 48
      Width = 49
      Height = 21
      TabOrder = 1
      Text = '0.001'
    end
    object AChangeEdit: TEdit
      Left = 136
      Top = 72
      Width = 49
      Height = 21
      TabOrder = 2
      Text = '0.1'
    end
    object WChangeEdit: TEdit
      Left = 136
      Top = 96
      Width = 49
      Height = 21
      TabOrder = 3
      Text = '0.1'
    end
  end
  object AutoRangeCheckBox: TCheckBox
    Left = 8
    Top = 127
    Width = 136
    Height = 17
    Caption = 'Set Range to Data'
    Checked = True
    State = cbChecked
    TabOrder = 3
  end
end
