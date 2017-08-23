object gStrainDialog: TgStrainDialog
  Left = 0
  Top = 0
  Caption = 'gStrain Dialog'
  ClientHeight = 91
  ClientWidth = 473
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object label1: TLabel
    Left = 40
    Top = 16
    Width = 147
    Height = 13
    Caption = 'Set the g-strain in percentage '
  end
  object Label2: TLabel
    Left = 280
    Top = 56
    Width = 11
    Height = 13
    Caption = '%'
  end
  object gStrainCheckBox: TCheckBox
    Left = 40
    Top = 56
    Width = 137
    Height = 17
    Caption = 'Use g-strain'
    TabOrder = 0
  end
  object gStrainEdit: TEdit
    Left = 200
    Top = 52
    Width = 65
    Height = 21
    TabOrder = 1
    Text = '0.0'
  end
  object BitBtn1: TBitBtn
    Left = 384
    Top = 48
    Width = 75
    Height = 25
    TabOrder = 2
    OnClick = BitBtn1Click
    Kind = bkOK
  end
end
