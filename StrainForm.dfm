object StrainDialog: TStrainDialog
  Left = 192
  Top = 107
  Width = 398
  Height = 223
  Caption = 'Strain Parameters'
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
    Width = 369
    Height = 145
    Caption = 'Crystal Field Strain Parameters'
    TabOrder = 0
    object Label1: TLabel
      Left = 208
      Top = 40
      Width = 21
      Height = 13
      Caption = 'GHz'
    end
    object Label2: TLabel
      Left = 208
      Top = 72
      Width = 21
      Height = 13
      Caption = 'GHz'
    end
    object B20StrainEdit: TEdit
      Left = 120
      Top = 32
      Width = 81
      Height = 21
      TabOrder = 0
      Text = '0.0'
    end
    object B20StrainCheckBox: TCheckBox
      Left = 16
      Top = 32
      Width = 81
      Height = 17
      Caption = 'B20 (D/3)'
      TabOrder = 1
    end
    object B22StrainCheckBox: TCheckBox
      Left = 16
      Top = 64
      Width = 57
      Height = 17
      Caption = 'B22 (E)'
      TabOrder = 2
    end
    object B22StrainEdit: TEdit
      Left = 120
      Top = 64
      Width = 81
      Height = 21
      TabOrder = 3
      Text = '0.0'
    end
  end
  object OKButton: TBitBtn
    Left = 304
    Top = 160
    Width = 75
    Height = 25
    TabOrder = 1
    OnClick = OKButtonClick
    Kind = bkOK
  end
  object CancelButton: TBitBtn
    Left = 216
    Top = 160
    Width = 75
    Height = 25
    TabOrder = 2
    Kind = bkCancel
  end
end
