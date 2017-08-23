object EPRSimulationDialog: TEPRSimulationDialog
  Left = 249
  Top = 154
  Hint = 'Number of points of simulation'
  AutoSize = True
  Caption = 'EPRSimulationDialog'
  ClientHeight = 449
  ClientWidth = 436
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBox4: TGroupBox
    Left = 380
    Top = 331
    Width = 55
    Height = 48
    Caption = 'Strain'
    TabOrder = 11
  end
  object SimParGroupBox: TGroupBox
    Left = 0
    Top = 0
    Width = 215
    Height = 210
    Caption = 'Simulation Parameters'
    TabOrder = 0
    object Label5: TLabel
      Left = 7
      Top = 26
      Width = 81
      Height = 13
      Caption = 'Number of Points'
    end
    object Label6: TLabel
      Left = 39
      Top = 52
      Width = 44
      Height = 13
      Caption = 'StartField'
    end
    object Label7: TLabel
      Left = 39
      Top = 78
      Width = 41
      Height = 13
      Caption = 'EndField'
    end
    object Label8: TLabel
      Left = 156
      Top = 52
      Width = 7
      Height = 13
      Caption = 'T'
    end
    object Label9: TLabel
      Left = 156
      Top = 78
      Width = 7
      Height = 13
      Caption = 'T'
    end
    object Label18: TLabel
      Left = 88
      Top = 136
      Width = 49
      Height = 13
      Caption = 'AngleStep'
    end
    object Label19: TLabel
      Left = 184
      Top = 128
      Width = 21
      Height = 13
      Caption = 'degr'
    end
    object Label22: TLabel
      Left = 32
      Top = 104
      Width = 56
      Height = 13
      Caption = 'Calcul order'
    end
    object Label10: TLabel
      Left = 104
      Top = 168
      Width = 28
      Height = 13
      Caption = 'Theta'
    end
    object Label12: TLabel
      Left = 104
      Top = 192
      Width = 15
      Height = 13
      Caption = 'Phi'
    end
    object Label29: TLabel
      Left = 184
      Top = 160
      Width = 21
      Height = 13
      Caption = 'degr'
    end
    object Label30: TLabel
      Left = 184
      Top = 184
      Width = 21
      Height = 13
      Caption = 'degr'
    end
    object Bevel1: TBevel
      Left = 16
      Top = 156
      Width = 185
      Height = 2
    end
    object NptsEdit: TEdit
      Left = 96
      Top = 16
      Width = 46
      Height = 21
      TabOrder = 0
      Text = '1024'
    end
    object StartFieldEdit: TEdit
      Left = 98
      Top = 45
      Width = 52
      Height = 21
      Hint = 'Lower Field limit (in Tesla)'
      TabOrder = 1
      Text = '8.5'
    end
    object EndFieldEdit: TEdit
      Left = 96
      Top = 72
      Width = 52
      Height = 21
      Hint = 'Upper Field limit (in Tesla)'
      TabOrder = 2
      Text = '8.62'
    end
    object AngleStepEdit: TEdit
      Left = 144
      Top = 128
      Width = 33
      Height = 21
      TabOrder = 3
      Text = '5'
      OnChange = AngleStepEditChange
    end
    object CalcOrderComboBox: TComboBox
      Left = 96
      Top = 96
      Width = 49
      Height = 21
      Style = csDropDownList
      ItemHeight = 13
      TabOrder = 4
      OnChange = CalcOrderComboBoxChange
      Items.Strings = (
        '1'
        '2'
        '3')
    end
    object PowderRadioButton: TRadioButton
      Left = 16
      Top = 136
      Width = 65
      Height = 17
      Caption = 'Powder'
      Checked = True
      TabOrder = 5
      TabStop = True
    end
    object CrystalRadioButton: TRadioButton
      Left = 16
      Top = 168
      Width = 57
      Height = 17
      Caption = 'Crystal'
      TabOrder = 6
    end
    object ThetaEdit: TEdit
      Left = 144
      Top = 160
      Width = 33
      Height = 21
      TabOrder = 7
      Text = '0'
    end
    object PhiEdit: TEdit
      Left = 144
      Top = 184
      Width = 33
      Height = 21
      TabOrder = 8
      Text = '0'
    end
  end
  object EPRLineGroupBox: TGroupBox
    Left = 221
    Top = 87
    Width = 215
    Height = 161
    Caption = 'Lineshape and width'
    TabOrder = 1
    object ModulationLabel: TLabel
      Left = 13
      Top = 114
      Width = 52
      Height = 13
      Caption = 'Modulation'
    end
    object Label1: TLabel
      Left = 65
      Top = 138
      Width = 15
      Height = 13
      Caption = 'mT'
    end
    object Label11: TLabel
      Left = 108
      Top = 15
      Width = 7
      Height = 13
      Caption = 'X'
    end
    object Label13: TLabel
      Left = 108
      Top = 67
      Width = 7
      Height = 13
      Caption = 'Z'
    end
    object Label14: TLabel
      Left = 108
      Top = 41
      Width = 7
      Height = 13
      Caption = 'Y'
    end
    object Phaselabel: TLabel
      Left = 88
      Top = 114
      Width = 54
      Height = 13
      Caption = 'Phase(deg)'
    end
    object wxLabel: TLabel
      Left = 192
      Top = 16
      Width = 15
      Height = 13
      Caption = 'mT'
    end
    object wyLabel: TLabel
      Left = 192
      Top = 40
      Width = 15
      Height = 13
      Caption = 'mT'
    end
    object wzLabel: TLabel
      Left = 192
      Top = 64
      Width = 15
      Height = 13
      Caption = 'mT'
    end
    object VoigtAlphaLabel: TLabel
      Left = 62
      Top = 95
      Width = 130
      Height = 13
      Caption = '% Lorentz (0:Gss...100:Lor.)'
    end
    object ModulationEdit: TEdit
      Left = 13
      Top = 132
      Width = 46
      Height = 21
      TabOrder = 0
      Text = '0.1'
    end
    object LwxEdit: TEdit
      Left = 121
      Top = 9
      Width = 65
      Height = 21
      TabOrder = 1
      Text = '1'
      OnChange = LwxEditChange
    end
    object LwyEdit: TEdit
      Left = 121
      Top = 35
      Width = 65
      Height = 21
      TabOrder = 2
      Text = '1'
      OnChange = LwyEditChange
    end
    object LwzEdit: TEdit
      Left = 121
      Top = 61
      Width = 65
      Height = 21
      TabOrder = 3
      Text = '1'
      OnChange = LwzEditChange
    end
    object DispersEdit: TEdit
      Left = 88
      Top = 132
      Width = 41
      Height = 21
      TabOrder = 4
      Text = '0'
      OnChange = DispersEditChange
    end
    object UpDown1: TUpDown
      Left = 129
      Top = 132
      Width = 15
      Height = 21
      Associate = DispersEdit
      Min = -180
      Max = 360
      TabOrder = 5
    end
    object DerivRadioButton: TRadioButton
      Left = 152
      Top = 120
      Width = 49
      Height = 17
      Caption = 'Deriv'
      Checked = True
      TabOrder = 6
      TabStop = True
    end
    object IntegrRadioButton: TRadioButton
      Left = 152
      Top = 136
      Width = 49
      Height = 17
      Caption = 'Integr'
      TabOrder = 7
    end
    object MixEdit: TEdit
      Left = 8
      Top = 88
      Width = 33
      Height = 21
      Enabled = False
      TabOrder = 8
      Text = '0'
    end
    object MixUpDown: TUpDown
      Left = 41
      Top = 88
      Width = 15
      Height = 21
      Associate = MixEdit
      TabOrder = 9
    end
    object ShapeButtonGroup: TRadioGroup
      Left = 3
      Top = 15
      Width = 99
      Height = 67
      Caption = 'Shape'
      ItemIndex = 1
      Items.Strings = (
        'Gauss'
        'Lorentz'
        '(pseudo)Voigt')
      TabOrder = 10
      OnClick = ShapeButtonGroupClick
    end
  end
  object GroupBox1: TGroupBox
    Left = 220
    Top = 0
    Width = 215
    Height = 81
    Caption = 'Experiment Parameters'
    TabOrder = 2
    object Label2: TLabel
      Left = 42
      Top = 26
      Width = 50
      Height = 13
      Caption = 'Frequency'
    end
    object Label3: TLabel
      Left = 29
      Top = 52
      Width = 60
      Height = 13
      Caption = 'Temperature'
    end
    object Label27: TLabel
      Left = 160
      Top = 24
      Width = 21
      Height = 13
      Caption = 'GHz'
    end
    object Label28: TLabel
      Left = 160
      Top = 48
      Width = 7
      Height = 13
      Caption = 'K'
    end
    object FreqEdit: TEdit
      Left = 101
      Top = 20
      Width = 52
      Height = 21
      TabOrder = 0
      Text = '240'
      OnChange = FreqEditChange
    end
    object TempEdit: TEdit
      Left = 101
      Top = 46
      Width = 52
      Height = 21
      TabOrder = 1
      Text = '290'
      OnChange = TempEditChange
    end
  end
  object BitBtn1: TBitBtn
    Left = 221
    Top = 385
    Width = 100
    Height = 27
    Caption = 'Load'
    TabOrder = 3
    Glyph.Data = {
      76010000424D7601000000000000760000002800000020000000100000000100
      04000000000000010000130B0000130B00001000000000000000000000000000
      800000800000008080008000000080008000808000007F7F7F00BFBFBF000000
      FF0000FF000000FFFF00FF000000FF00FF00FFFF0000FFFFFF0033333333B333
      333B33FF33337F3333F73BB3777BB7777BB3377FFFF77FFFF77333B000000000
      0B3333777777777777333330FFFFFFFF07333337F33333337F333330FFFFFFFF
      07333337F3FF3FFF7F333330F00F000F07333337F77377737F333330FFFFFFFF
      07333FF7F3FFFF3F7FFFBBB0F0000F0F0BB37777F7777373777F3BB0FFFFFFFF
      0BBB3777F3FF3FFF77773330F00F000003333337F773777773333330FFFF0FF0
      33333337F3FF7F37F3333330F08F0F0B33333337F7737F77FF333330FFFF003B
      B3333337FFFF77377FF333B000000333BB33337777777F3377FF3BB3333BB333
      3BB33773333773333773B333333B3333333B7333333733333337}
    NumGlyphs = 2
  end
  object BitBtn2: TBitBtn
    Left = 327
    Top = 385
    Width = 100
    Height = 27
    Caption = 'Save'
    TabOrder = 4
    Glyph.Data = {
      76010000424D7601000000000000760000002800000020000000100000000100
      04000000000000010000130B0000130B00001000000000000000000000000000
      800000800000008080008000000080008000808000007F7F7F00BFBFBF000000
      FF0000FF000000FFFF00FF000000FF00FF00FFFF0000FFFFFF00333333330070
      7700333333337777777733333333008088003333333377F73377333333330088
      88003333333377FFFF7733333333000000003FFFFFFF77777777000000000000
      000077777777777777770FFFFFFF0FFFFFF07F3333337F3333370FFFFFFF0FFF
      FFF07F3FF3FF7FFFFFF70F00F0080CCC9CC07F773773777777770FFFFFFFF039
      99337F3FFFF3F7F777F30F0000F0F09999937F7777373777777F0FFFFFFFF999
      99997F3FF3FFF77777770F00F000003999337F773777773777F30FFFF0FF0339
      99337F3FF7F3733777F30F08F0F0337999337F7737F73F7777330FFFF0039999
      93337FFFF7737777733300000033333333337777773333333333}
    NumGlyphs = 2
  end
  object BitBtn4: TBitBtn
    Left = 327
    Top = 416
    Width = 100
    Height = 27
    TabOrder = 5
    OnClick = BitBtn4Click
    Kind = bkOK
  end
  object BitBtn5: TBitBtn
    Left = 221
    Top = 416
    Width = 100
    Height = 27
    TabOrder = 6
    OnClick = BitBtn5Click
    Kind = bkCancel
  end
  object PopGroupBox: TGroupBox
    Left = 0
    Top = 210
    Width = 215
    Height = 135
    Caption = 'Population'
    TabOrder = 7
    object Label20: TLabel
      Left = 88
      Top = 32
      Width = 71
      Height = 13
      Caption = '(Px-Pz)/(Py-Pz)'
    end
    object Label23: TLabel
      Left = 16
      Top = 74
      Width = 29
      Height = 13
      Caption = 'Px/Pz'
    end
    object Label24: TLabel
      Left = 112
      Top = 72
      Width = 29
      Height = 13
      Caption = 'Py/Pz'
    end
    object BoltzmannCheckBox: TCheckBox
      Left = 8
      Top = 15
      Width = 78
      Height = 13
      Caption = 'Boltzmann'
      TabOrder = 0
    end
    object GroupBox3: TGroupBox
      Left = 8
      Top = 89
      Width = 201
      Height = 40
      Caption = 'Triplet pops'
      TabOrder = 1
      object Label15: TLabel
        Left = 3
        Top = 16
        Width = 12
        Height = 13
        Caption = 'Px'
      end
      object Label16: TLabel
        Left = 67
        Top = 16
        Width = 12
        Height = 13
        Caption = 'Py'
      end
      object Label17: TLabel
        Left = 131
        Top = 16
        Width = 12
        Height = 13
        Caption = 'Pz'
      end
      object PxEdit: TEdit
        Left = 20
        Top = 12
        Width = 46
        Height = 21
        Enabled = False
        TabOrder = 0
        Text = '0'
        OnChange = PxEditChange
      end
      object PzEdit: TEdit
        Left = 148
        Top = 12
        Width = 46
        Height = 21
        TabOrder = 1
        Text = '1'
      end
      object PyEdit: TEdit
        Left = 84
        Top = 12
        Width = 46
        Height = 21
        Enabled = False
        TabOrder = 2
        Text = '0'
      end
    end
    object PopRatioEdit: TEdit
      Left = 160
      Top = 24
      Width = 49
      Height = 21
      TabOrder = 2
      Text = '1'
      OnChange = PopRatioEditChange
    end
    object XZratioEdit: TEdit
      Left = 47
      Top = 67
      Width = 58
      Height = 21
      TabOrder = 3
      Text = '1'
      OnChange = XZratioEdit1Change
    end
    object SOISCcheckbox: TCheckBox
      Left = 8
      Top = 30
      Width = 57
      Height = 17
      Caption = 'SO-ISC'
      TabOrder = 4
    end
    object RPISCcheckbox: TCheckBox
      Left = 8
      Top = 45
      Width = 73
      Height = 17
      Caption = 'RP-ISC'
      TabOrder = 5
    end
    object YZratioEdit: TEdit
      Left = 152
      Top = 64
      Width = 57
      Height = 21
      TabOrder = 6
      Text = '1'
      OnChange = YZratioEditChange
    end
  end
  object StrainEdit: TEdit
    Left = 386
    Top = 353
    Width = 41
    Height = 21
    TabOrder = 8
    Text = '0'
  end
  object PartialOrientationGroupBox: TGroupBox
    Left = 220
    Top = 248
    Width = 215
    Height = 81
    Caption = 'Partial Orientation'
    TabOrder = 9
    object Label33: TLabel
      Left = 16
      Top = 56
      Width = 14
      Height = 16
      Caption = 'J '
      Font.Charset = SYMBOL_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Symbol'
      Font.Style = [fsBold]
      ParentFont = False
    end
    object Label34: TLabel
      Left = 88
      Top = 56
      Width = 6
      Height = 16
      Caption = 'f'
      Font.Charset = SYMBOL_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Symbol'
      Font.Style = [fsBold]
      ParentFont = False
    end
    object Non_OrientedRadioButton: TRadioButton
      Left = 8
      Top = 16
      Width = 113
      Height = 17
      Caption = 'No orientation'
      Checked = True
      TabOrder = 0
      TabStop = True
    end
    object OrientedRadioButton: TRadioButton
      Left = 8
      Top = 32
      Width = 145
      Height = 17
      Caption = 'Oriented, Energy (in kT)'
      TabOrder = 1
    end
    object PartOrThetaEdit: TEdit
      Left = 32
      Top = 56
      Width = 41
      Height = 21
      TabOrder = 2
      Text = '0.0'
    end
    object PartOrPhiEdit: TEdit
      Left = 104
      Top = 56
      Width = 41
      Height = 21
      TabOrder = 3
      Text = '0.0'
    end
    object OrderParEdit: TEdit
      Left = 144
      Top = 24
      Width = 41
      Height = 21
      TabOrder = 4
      Text = '1'
    end
  end
  object GroupBox2: TGroupBox
    Left = 0
    Top = 344
    Width = 215
    Height = 105
    Caption = 'Decay and SLR Rates (10^6 s-1)'
    TabOrder = 10
    object Label35: TLabel
      Left = 8
      Top = 56
      Width = 31
      Height = 13
      Caption = 'Decay'
    end
    object Label36: TLabel
      Left = 96
      Top = 32
      Width = 7
      Height = 13
      Caption = 'X'
    end
    object Label37: TLabel
      Left = 144
      Top = 32
      Width = 7
      Height = 13
      Caption = 'Y'
    end
    object Label38: TLabel
      Left = 192
      Top = 32
      Width = 7
      Height = 13
      Caption = 'Z'
    end
    object Label39: TLabel
      Left = 8
      Top = 80
      Width = 43
      Height = 13
      Caption = 'SL Relax'
    end
    object Label40: TLabel
      Left = 8
      Top = 16
      Width = 58
      Height = 13
      Caption = 'Time (10-6s)'
    end
    object kxEdit: TEdit
      Left = 72
      Top = 48
      Width = 40
      Height = 21
      TabOrder = 0
      Text = '0'
    end
    object kyEdit: TEdit
      Left = 120
      Top = 48
      Width = 40
      Height = 21
      TabOrder = 1
      Text = '0'
    end
    object kzEdit: TEdit
      Left = 168
      Top = 48
      Width = 40
      Height = 21
      TabOrder = 2
      Text = '0'
    end
    object SLRxEdit: TEdit
      Left = 72
      Top = 72
      Width = 40
      Height = 21
      TabOrder = 3
      Text = '0'
    end
    object SLRyEdit: TEdit
      Left = 120
      Top = 72
      Width = 40
      Height = 21
      TabOrder = 4
      Text = '0'
    end
    object SLRzEdit: TEdit
      Left = 168
      Top = 72
      Width = 40
      Height = 21
      TabOrder = 5
      Text = '0'
    end
    object TimeEdit: TEdit
      Left = 8
      Top = 32
      Width = 49
      Height = 21
      TabOrder = 6
      Text = '0'
    end
    object FixedbCheckBox: TCheckBox
      Left = 152
      Top = 16
      Width = 97
      Height = 17
      Caption = 'Fixed b'
      TabOrder = 7
    end
  end
  object PolarizationGroupBox: TGroupBox
    Left = 221
    Top = 331
    Width = 153
    Height = 48
    Caption = 'mm-wave Polarization'
    TabOrder = 12
    object PolarizationListBox: TComboBox
      Left = 5
      Top = 22
      Width = 145
      Height = 21
      Style = csDropDownList
      ItemHeight = 13
      ItemIndex = 0
      TabOrder = 0
      Text = 'in xy plane, perp B0'
      Items.Strings = (
        'in xy plane, perp B0'
        'in B0-Z plane, perp B0'
        'along B0'
        'along x'
        'along y'
        'along z'
        'right circular'
        'left circular')
    end
  end
end
