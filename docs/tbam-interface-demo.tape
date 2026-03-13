Output docs/tbam-interface-demo.gif

Set Shell "zsh"
Set Width 1800
Set Height 900
Set FontSize 16
Set TypingSpeed 80ms
Set Theme "Builtin Solarized Dark"

Type "cd /Users/parsaghadermarzi/Desktop/Academics/Projects/TBAM"
Enter
Sleep 500ms

Type "source .venv/bin/activate"
Enter
Sleep 500ms

Type "python -m tbam.cli demo-data --out-dir sample_data --reads 450"
Enter
Sleep 2s

Type "python -m tbam.cli view sample_data/demo.bam --reference sample_data/demo.fa --window-bp 40 --panel-reads 35"
Enter
Sleep 2s

Type "d"
Sleep 500ms
Type "d"
Sleep 500ms
Type "]"
Sleep 500ms
Type "a"
Sleep 500ms
Type "0"
Sleep 900ms
Type "q"
Sleep 700ms
