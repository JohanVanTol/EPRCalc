object EnergyForm: TEnergyForm
  Left = 0
  Top = 0
  Caption = 'Energy Level Calculation'
  ClientHeight = 286
  ClientWidth = 426
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
    Left = 52
    Top = 67
    Width = 22
    Height = 13
    Caption = 'Field'
  end
  object Label2: TLabel
    Left = 167
    Top = 67
    Width = 6
    Height = 13
    Caption = 'T'
  end
  object FixedFieldRadioButton: TRadioButton
    Left = 40
    Top = 17
    Width = 161
    Height = 17
    Caption = 'Orientation dependent'
    TabOrder = 0
  end
  object FixedOrientRadioButton: TRadioButton
    Left = 248
    Top = 17
    Width = 113
    Height = 17
    Caption = 'Field dependent'
    TabOrder = 1
  end
  object Edit1: TEdit
    Left = 80
    Top = 64
    Width = 81
    Height = 21
    TabOrder = 2
    Text = 'Edit1'
  end
  object Edit2: TEdit
    Left = 280
    Top = 64
    Width = 81
    Height = 21
    TabOrder = 3
    Text = 'Edit1'
  end
end
