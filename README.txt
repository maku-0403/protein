***パスをコピーするときは全て/になっているか確認する。\だと動作しない。***

1. SSDのパス名をコピー

2. ターミナル（黒い画面に文字しかないやつ）で
cd コピーしたパス
を入力してEnterを押す。

3. 次に
rsync -avz --delete rsync.pdbj.org::ftp_data/structures/divided/mmCIF/ ./mmCIF
を入力してEnterを押す。

4. 全て処理が終わったらSSD内に新しくできたmmCIFのファイルのパス名をコピーする

5. ターミナルで
cd コピーしたパス
でEnter、続けて
find . -type f -name "*.gz" -exec gunzip {} +
でEnterして、終了まで待つ（結構時間かかる）

6. protein -> for_cif -> cif_to_csv_complete_cosw.pyを開く。

7. SSDに"CSV"というファイルを作る。

8. root_dirにSSD内のmmCIFファイルのパス、out_dirに7で作ったCSVファイルのパスを貼り付ける。

9. cif_to_csv_complete_cosw.pyを実行して終了するまで待つ（結構時間かかる）

10. SSDに"sphere"というファイルを作る。

11. root_dirにSSD内のmmCIFファイルのパス、out_dirに10で作ったsphereファイルのパスを貼り付ける。

12. cif_temperature_sphere.pyを実行して終了するまで待つ（結構時間かかる）

13. SSDに"upper_25"というファイルを作る。

14. root_dirにSSD内のsphereファイルのパス、out_dirに13で作ったupper_25ファイルのパスを貼り付ける。

15. csv_retouch_upper.py内で切る割合を適宜設定する。

16. csv_retouch_upper.pyを実行して終了するまで待つ（結構時間かかる）

***パスをコピーするときは全て/になっているか確認する。\だと動作しない。***

find ./mmCIF -type f -name "*.gz" -exec gunzip {} +
