object gStrainDialog: TgStrainDialog
  Left = 0
  Top = 0
  Caption = 'g Strain Dialog'
  ClientHeight = 75
  ClientWidth = 373
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Label1: TLabel
    Left = 207
    Top = 43
    Width = 17
    Height = 13
    Caption = ' % '
  end
  object Label2: TLabel
    Left = 40
    Top = 8
    Width = 100
    Height = 13
    Caption = 'Set the g-strain in %'
  end
  object gStrainCheckBox: TCheckBox
    Left = 40
    Top = 40
    Width = 97
    Height = 17
    Caption = 'Use g Strain'
    TabOrder = 0
  end
  object gStrainEdit: TEdit
    Left = 152
    Top = 40
    Width = 49
    Height = 21
    TabOrder = 1
    Text = '0.0'
  end
  object BitBtn1: TBitBtn
    Left = 288
    Top = 36
    Width = 75
    Height = 25
    TabOrder = 2
    OnClick = BitBtn1Click
    Kind = bkOK
  end
end
