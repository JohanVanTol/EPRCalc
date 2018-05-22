object ExponentialDecayForm: TExponentialDecayForm
  Left = 0
  Top = 0
  Caption = 'ExponentialFit'
  ClientHeight = 442
  ClientWidth = 556
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object OKButton: TBitBtn
    Left = 464
    Top = 409
    Width = 75
    Height = 25
    Caption = 'Done'
    TabOrder = 0
    OnClick = OKButtonClick
    Kind = bkOK
  end
  object CancelButton: TBitBtn
    Left = 383
    Top = 409
    Width = 75
    Height = 25
    TabOrder = 1
    OnClick = CancelButtonClick
    Kind = bkCancel
  end
  object SimulateButton: TButton
    Left = 23
    Top = 409
    Width = 75
    Height = 25
    Caption = 'Simulate'
    TabOrder = 2
    OnClick = SimulateButtonClick
  end
  object StartTimeEdit: TEdit
    Left = 104
    Top = 16
    Width = 121
    Height = 21
    TabOrder = 3
    Text = 'StartTimeEdit'
  end
  object EndTimeEdit: TEdit
    Left = 240
    Top = 16
    Width = 121
    Height = 21
    TabOrder = 4
    Text = 'EndTimeEdit'
  end
  object SingleExpGroupBox: TGroupBox
    Left = 23
    Top = 47
    Width = 515
    Height = 346
    Caption = 'Exponential Decay'
    TabOrder = 5
    object y0CheckBox: TCheckBox
      Left = 32
      Top = 32
      Width = 97
      Height = 17
      Caption = 'Baseline'
      TabOrder = 0
    end
    object Amp1CheckBox: TCheckBox
      Left = 32
      Top = 71
      Width = 97
      Height = 17
      Caption = 'Amplitude 1'
      TabOrder = 1
    end
    object Exp1CheckBox: TCheckBox
      Left = 32
      Top = 94
      Width = 97
      Height = 17
      Caption = 'Decay Time 1'
      TabOrder = 2
    end
    object Amp2ErrorEdit: TEdit
      Left = 352
      Top = 129
      Width = 121
      Height = 21
      Enabled = False
      TabOrder = 3
      Text = '0'
    end
    object Amp2CheckBox: TCheckBox
      Left = 32
      Top = 131
      Width = 97
      Height = 17
      Caption = 'Amplitude 2'
      Enabled = False
      TabOrder = 4
    end
    object Exp2CheckBox: TCheckBox
      Left = 32
      Top = 152
      Width = 97
      Height = 17
      Caption = 'Decay Time 2'
      Enabled = False
      TabOrder = 5
    end
    object GaussAmpCheckBox: TCheckBox
      Left = 32
      Top = 194
      Width = 97
      Height = 17
      Caption = 'Gauss Amplitude'
      Enabled = False
      TabOrder = 6
    end
    object GaussZeroCheckBox: TCheckBox
      Left = 32
      Top = 217
      Width = 97
      Height = 17
      Caption = 'Time zero'
      Enabled = False
      TabOrder = 7
    end
    object GaussSigmaCheckBox: TCheckBox
      Left = 32
      Top = 240
      Width = 123
      Height = 17
      Caption = 'Gaussian HalfWidth'
      Enabled = False
      TabOrder = 8
    end
    object Amp2Edit: TEdit
      Left = 176
      Top = 129
      Width = 121
      Height = 21
      Enabled = False
      TabOrder = 9
      Text = '0'
    end
    object Exp2Edit: TEdit
      Left = 176
      Top = 150
      Width = 121
      Height = 21
      Enabled = False
      TabOrder = 10
      Text = '1'
    end
    object GaussAmpEdit: TEdit
      Left = 176
      Top = 190
      Width = 121
      Height = 21
      Enabled = False
      TabOrder = 11
      Text = '0'
    end
    object GaussZeroEdit: TEdit
      Left = 176
      Top = 215
      Width = 121
      Height = 21
      Enabled = False
      TabOrder = 12
      Text = '0'
    end
    object GaussSigmaEdit: TEdit
      Left = 176
      Top = 236
      Width = 121
      Height = 21
      Enabled = False
      TabOrder = 13
      Text = '1'
    end
    object Exp2ErrorEdit: TEdit
      Left = 352
      Top = 150
      Width = 121
      Height = 21
      Enabled = False
      TabOrder = 14
      Text = '0'
    end
    object GaussAmpErrorEdit: TEdit
      Left = 352
      Top = 190
      Width = 121
      Height = 21
      Enabled = False
      TabOrder = 15
      Text = '0'
    end
    object GaussZeroErrorEdit: TEdit
      Left = 352
      Top = 213
      Width = 121
      Height = 21
      Enabled = False
      TabOrder = 16
      Text = '0'
    end
    object GaussSigmaErrorEdit: TEdit
      Left = 352
      Top = 236
      Width = 121
      Height = 21
      Enabled = False
      TabOrder = 17
      Text = '0'
    end
    object BiRadioButton: TRadioButton
      Left = 16
      Top = 326
      Width = 113
      Height = 17
      Caption = 'Bi-exponential'
      TabOrder = 18
      OnClick = BiRadioButtonClick
    end
    object MonoGaussRadioButton: TRadioButton
      Left = 184
      Top = 324
      Width = 137
      Height = 17
      Caption = 'Mono-Exp + Gaussian'
      TabOrder = 19
      OnClick = MonoGaussRadioButtonClick
    end
    object MonoRadioButton: TRadioButton
      Left = 16
      Top = 303
      Width = 113
      Height = 17
      Caption = 'Mono-exponential'
      Checked = True
      TabOrder = 20
      TabStop = True
      OnClick = MonoRadioButtonClick
    end
    object StretchedRadioButton: TRadioButton
      Left = 184
      Top = 301
      Width = 113
      Height = 17
      Caption = 'Stretched Exponential'
      TabOrder = 21
      OnClick = StretchedRadioButtonClick
    end
    object StretchCheckBox: TCheckBox
      Left = 32
      Top = 280
      Width = 113
      Height = 17
      Caption = 'Stretch Parameter'
      Enabled = False
      TabOrder = 22
    end
    object StretchEdit: TEdit
      Left = 176
      Top = 274
      Width = 121
      Height = 21
      Enabled = False
      TabOrder = 23
      Text = '1'
    end
    object StretchErrorEdit: TEdit
      Left = 352
      Top = 274
      Width = 121
      Height = 21
      Enabled = False
      TabOrder = 24
      Text = '0'
    end
  end
  object y0Edit: TEdit
    Left = 199
    Top = 77
    Width = 121
    Height = 21
    TabOrder = 6
    Text = 'y0Edit'
  end
  object Amp1Edit: TEdit
    Left = 199
    Top = 116
    Width = 121
    Height = 21
    TabOrder = 7
    Text = 'Amp1Edit'
  end
  object Exp1Edit: TEdit
    Left = 199
    Top = 139
    Width = 121
    Height = 21
    TabOrder = 8
    Text = 'Exp1Edit'
  end
  object y0ErrorEdit: TEdit
    Left = 375
    Top = 77
    Width = 121
    Height = 21
    TabOrder = 9
    Text = 'y0ErrorEdit'
  end
  object Amp1ErrorEdit: TEdit
    Left = 375
    Top = 116
    Width = 121
    Height = 21
    TabOrder = 10
    Text = 'Amp1ErrorEdit'
  end
  object Exp1ErrorEdit: TEdit
    Left = 375
    Top = 139
    Width = 121
    Height = 21
    TabOrder = 11
    Text = 'Exp1ErrorEdit'
  end
  object NptsEdit: TEdit
    Left = 32
    Top = 16
    Width = 51
    Height = 21
    TabOrder = 12
    Text = 'NptsEdit'
  end
  object CorrectBaseLineButton: TButton
    Left = 440
    Top = 16
    Width = 99
    Height = 25
    Caption = 'Correct BaseLine'
    TabOrder = 13
    OnClick = CorrectBaseLineButtonClick
  end
  object FitButton: TButton
    Left = 104
    Top = 409
    Width = 75
    Height = 25
    Caption = 'Fit 1 cycle'
    TabOrder = 14
    OnClick = Fit10ButtonClick
  end
  object Button1: TButton
    Left = 185
    Top = 409
    Width = 75
    Height = 25
    Caption = 'Fit 10 cycles'
    TabOrder = 15
    OnClick = Fit10ButtonClick
  end
end
