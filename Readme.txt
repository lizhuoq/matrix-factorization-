功能：1.求方阵A的特征值
2.对方阵A进行LU分解，对矩阵A进行QR（Gram-Schmidt）、Householder和Givens分解，并在此基础上求Ax=b的解
3.对矩阵A进行URV分解
语言：python
依赖库：numpy
1.factoriza.py
pp_lu函数：采用的方法为部分主元LU分解，其中输入为A矩阵（要求A为方阵），xrow为某一列绝对值最大元素所在行，输出为P为行变换矩阵，L为对角线为1的下三角矩阵，U为上三角矩阵。
qr_gs函数：gram-schmidt法将A矩阵变为单位正交矩阵Q左乘上三角阵R，其中r为A矩阵某一列向量在之前已经确定好的方向向量上的投影，然后该列向量减去在已经确定好的方向上的投影向量，anorm表示新的正交方向的向量的模长，q为单位正交向量，输出Q为单位正交矩阵，R为上三角矩阵。
householder函数：输入A为一个矩阵，verbose控制是否打印P矩阵或Q矩阵，a为A矩阵位于对角线下方的列向量，e为和a长度相同，首个元素为a的模长，其余元素为0的向量，输出P矩阵和T矩阵使得PA=T，当A为方阵时，T为上三角矩阵，P为正交矩阵。
givens函数：输入A为一个矩阵，r和c为矩阵主对角线下部分的元素位置，s为主对角线上的元素与其下行中非零元素的平方和，主对角线上的元素除s为旋转角度的cos值，其下非零元素除s为旋转角度的sin值，输出P为正交矩阵，使得PA=T，当A为方阵时，T为上三角矩阵。
urv函数：输入A为一个矩阵，P为A做householder reduction的一个正交矩阵，B为使得PA=B，取B矩阵的前非零行，Q为做householder reduction的一个正交矩阵，使得QB.T=T，Q为一正交矩阵，取T的前非零行作为T，最后输出U为P.T，V为Q.T，R为和A矩阵相同大小，左上角为T.T，其余元素为0的矩阵。
factorization函数：输入A矩阵，输入矩阵分解方法，输出为对应的分解矩阵。
2.solution.py
sub_matrix函数：输入为A为一个方阵，r和c分别为行和列的位置，verbose控制是否打印子矩阵，输出SM矩阵为去掉r行和c列的一个子矩阵。
det函数：输入A为一个方阵，verbose控制是否打印行列式值，对于第一行每个元素求代数余子式，并对乘积求和，输出det_value为行列式值
ly_pb函数：输入L为一个下三角阵，b为一个长度等于L行数的列向量，P默认值为单位阵，当作PA=LU分解时，可以输入P矩阵，按行从上到下依次消去，最后输出为列向量y。
ux_y函数：输入U 为一个下三角阵，y为一个长度等于U行数的列向量，按行从下到上依次消去，最后输出为列向量x。
qr_resolution函数：输入Q和R为经过gram-Schmidt分解得到的矩阵，b为等于Q的行数的列向量，把R的前非零行赋值给R，c为Q.T左乘b，令c的元素个数等于R的行数，R为一个上三角阵，c为等于R行数的列向量，把R和c输入ux_y函数，得到x为Ax=b的解。
householder_givens_solution函数：输入为经过householder reduction或givens reduction后的P和T矩阵，b为一列向量，令R为T前非零行矩阵，c为P左乘b，令c元素个数等于R的行数，R为一上三角阵，把R，c输入ux_y函数，得到x为Ax=b的解。
solution函数：输入A为一个矩阵，method为矩阵分解方法，b为一个列向量，如果A为方阵，打印特征值，输出相应的分解矩阵和Ax=b的解向量x。



