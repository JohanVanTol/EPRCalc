object FitCycleDialog: TFitCycleDialog
  Left = 192
  Top = 107
  AutoSize = True
  Caption = 'Fit cycles'
  ClientHeight = 313
  ClientWidth = 441
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object HamiltonianFitGroupBox: TGroupBox
    Left = 0
    Top = 1
    Width = 140
    Height = 225
    Caption = 'ZFS and g'
    TabOrder = 0
    object DFitCheckBox: TCheckBox
      Left = 16
      Top = 16
      Width = 33
      Height = 17
      Caption = 'D'
      TabOrder = 0
    end
    object EFitCheckBox: TCheckBox
      Left = 16
      Top = 32
      Width = 33
      Height = 17
      Caption = 'E'
      TabOrder = 1
    end
    object GxxFitCheckBox: TCheckBox
      Left = 16
      Top = 80
      Width = 50
      Height = 17
      Caption = 'Gxx'
      TabOrder = 2
    end
    object GyyFitCheckBox: TCheckBox
      Left = 16
      Top = 104
      Width = 48
      Height = 17
      Caption = 'Gyy'
      TabOrder = 3
    end
    object GzzFitCheckBox: TCheckBox
      Left = 16
      Top = 128
      Width = 48
      Height = 17
      Caption = 'Gzz'
      TabOrder = 4
    end
    object GxyFitCheckBox: TCheckBox
      Left = 16
      Top = 152
      Width = 48
      Height = 17
      Caption = 'Gxy'
      TabOrder = 5
    end
    object GxzFitCheckBox: TCheckBox
      Left = 16
      Top = 176
      Width = 48
      Height = 17
      Caption = 'Gxz'
      TabOrder = 6
    end
    object GyzFitCheckBox: TCheckBox
      Left = 16
      Top = 200
      Width = 48
      Height = 17
      Caption = 'Gyz'
      TabOrder = 7
    end
    object GisoCheckBox: TCheckBox
      Left = 16
      Top = 56
      Width = 50
      Height = 17
      Caption = 'Giso'
      TabOrder = 8
    end
    object DStrainFitCheckBox: TCheckBox
      Left = 48
      Top = 16
      Width = 33
      Height = 17
      Caption = 'St'
      TabOrder = 9
    end
    object EStrainFitCheckBox: TCheckBox
      Left = 48
      Top = 32
      Width = 33
      Height = 17
      Caption = 'St'
      TabOrder = 10
    end
    object B40FitCheckBox: TCheckBox
      Left = 72
      Top = 55
      Width = 45
      Height = 17
      Caption = 'B40'
      TabOrder = 11
    end
  end
  object GroupBox2: TGroupBox
    Left = 233
    Top = 0
    Width = 121
    Height = 113
    Caption = 'Populations'
    TabOrder = 1
    object AbsDisFitCheckBox: TCheckBox
      Left = 16
      Top = 16
      Width = 97
      Height = 17
      Caption = 'Phase'
      TabOrder = 0
    end
    object ISCPopRatioFitCheckBox: TCheckBox
      Left = 16
      Top = 40
      Width = 89
      Height = 17
      Caption = 'ISC Pop. Ratio'
      TabOrder = 1
    end
    object XZratioFitCheckBox: TCheckBox
      Left = 16
      Top = 65
      Width = 97
      Height = 17
      Caption = 'X/Z Pop. Ratio'
      TabOrder = 2
    end
    object YZratioFitCheckBox: TCheckBox
      Left = 16
      Top = 88
      Width = 97
      Height = 17
      Caption = 'Y/Z Pop. Ratio'
      TabOrder = 3
    end
  end
  object GroupBox3: TGroupBox
    Left = 233
    Top = 113
    Width = 121
    Height = 113
    Caption = 'Linewidth'
    TabOrder = 2
    object WxFitCheckBox: TCheckBox
      Left = 16
      Top = 16
      Width = 97
      Height = 17
      Caption = 'Wx'
      TabOrder = 0
    end
    object WyFitCheckBox: TCheckBox
      Left = 16
      Top = 40
      Width = 97
      Height = 17
      Caption = 'Wy'
      TabOrder = 1
    end
    object WzFitCheckBox: TCheckBox
      Left = 16
      Top = 64
      Width = 97
      Height = 17
      Caption = 'Wz'
      TabOrder = 2
    end
    object WisoFitCheckBox: TCheckBox
      Left = 16
      Top = 90
      Width = 97
      Height = 17
      Caption = 'Wiso'
      TabOrder = 3
    end
  end
  object GroupBox1: TGroupBox
    Left = 0
    Top = 232
    Width = 185
    Height = 81
    Caption = 'FitParam'
    TabOrder = 3
    object Ncycles: TLabel
      Left = 24
      Top = 24
      Width = 38
      Height = 13
      Caption = 'Ncycles'
    end
    object Label56: TLabel
      Left = 160
      Top = 48
      Width = 8
      Height = 13
      Caption = '%'
    end
    object Label1: TLabel
      Left = 24
      Top = 56
      Width = 81
      Height = 13
      Caption = 'Min Improvement'
    end
    object NCycleEdit: TEdit
      Left = 112
      Top = 16
      Width = 41
      Height = 21
      TabOrder = 0
      Text = '1'
    end
    object ImprovementEdit: TEdit
      Left = 112
      Top = 48
      Width = 41
      Height = 21
      TabOrder = 1
      Text = '1'
    end
    object UpDown1: TUpDown
      Left = 153
      Top = 16
      Width = 15
      Height = 21
      Associate = NCycleEdit
      Min = 1
      Position = 1
      TabOrder = 2
    end
  end
  object BitBtn1: TBitBtn
    Left = 344
    Top = 241
    Width = 97
    Height = 33
    TabOrder = 4
    OnClick = BitBtn1Click
    Kind = bkCancel
  end
  object BitBtn2: TBitBtn
    Left = 344
    Top = 280
    Width = 97
    Height = 33
    TabOrder = 5
    OnClick = BitBtn2Click
    Kind = bkOK
  end
  object GroupBox4: TGroupBox
    Left = 146
    Top = 1
    Width = 81
    Height = 225
    Caption = 'Hyperfine'
    TabOrder = 6
    object AxxFitCheckBox: TCheckBox
      Left = 8
      Top = 80
      Width = 41
      Height = 17
      Caption = 'Axx'
      TabOrder = 0
    end
    object AyyFitCheckBox: TCheckBox
      Left = 8
      Top = 104
      Width = 49
      Height = 17
      Caption = 'Ayy'
      TabOrder = 1
    end
    object AzzFitCheckBox: TCheckBox
      Left = 8
      Top = 128
      Width = 49
      Height = 17
      Caption = 'Azz'
      TabOrder = 2
    end
    object AxyFitCheckBox: TCheckBox
      Left = 8
      Top = 152
      Width = 49
      Height = 17
      Caption = 'Axy'
      TabOrder = 3
    end
    object AxzFitCheckBox: TCheckBox
      Left = 8
      Top = 176
      Width = 49
      Height = 17
      Caption = 'Axy'
      TabOrder = 4
    end
    object AyzFitCheckBox: TCheckBox
      Left = 8
      Top = 200
      Width = 49
      Height = 17
      Caption = 'Ayz'
      TabOrder = 5
    end
    object AisoFitCheckBox: TCheckBox
      Left = 8
      Top = 56
      Width = 41
      Height = 17
      Caption = 'Aiso'
      TabOrder = 6
    end
  end
  object GroupBox5: TGroupBox
    Left = 360
    Top = 0
    Width = 73
    Height = 225
    Caption = 'Decay'
    TabOrder = 7
    object PzFitCheckBox: TCheckBox
      Left = 8
      Top = 16
      Width = 40
      Height = 17
      Caption = 'Pz'
      TabOrder = 0
    end
    object KxFitCheckBox: TCheckBox
      Left = 8
      Top = 64
      Width = 40
      Height = 17
      Caption = 'kx'
      TabOrder = 1
    end
    object KyFitCheckBox: TCheckBox
      Left = 8
      Top = 88
      Width = 40
      Height = 17
      Caption = 'ky'
      TabOrder = 2
    end
    object KzFitCheckBox: TCheckBox
      Left = 8
      Top = 112
      Width = 40
      Height = 17
      Caption = 'kz'
      TabOrder = 3
    end
    object SLRxFitCheckBox: TCheckBox
      Left = 8
      Top = 144
      Width = 50
      Height = 17
      Caption = 'SLR x'
      TabOrder = 4
    end
    object SLRyFitCheckBox: TCheckBox
      Left = 8
      Top = 168
      Width = 50
      Height = 17
      Caption = 'SLR y'
      TabOrder = 5
    end
    object SLRzFitCheckBox: TCheckBox
      Left = 8
      Top = 192
      Width = 50
      Height = 17
      Caption = 'SLR z'
      TabOrder = 6
    end
  end
end
