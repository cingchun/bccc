使用手册(for windows)

首先，确定你电脑有curl（现在大部分电脑都已自带curl）。
    如何确定? 
	    打开windows命令工具 cmd，输入 curl --help
	        若是curl各项参数说明表明你电脑有curl，若报错表明没有
    没有curl请自行安装 https://www.cnblogs.com/gered/p/10682298.html
第二，修改。
    记事本打开login.bat，把其中的username和password填入自己的学号和密码
最后，配置。
    将两个脚本加入windows计划任务


有很多同学咨询windows添加计划任务
    可参考 https://jingyan.baidu.com/article/154b463130041128ca8f41c7.html    
可能有些同学觉得点来点去设置比较麻烦，这里还有另一种方式
    直接将两个脚本放到这个文件夹下 C:\Users\admin\AppData\Roaming\Microsoft\Windows\Start Menu\Programs\Startup
    重启电脑 
    
姑且称这种为方式二，上一种为方式一
两种方式各有优劣：
方式一是永久加入，等回学校后，如果你不想使用，右键禁用即可；若以后再想用，比如出、差放暑假，右键启用。
    这样想用或不用比较方便
方式二临时使用，如不想使用要去文件夹将其删除，再想使用时，再从其它地方拷贝到文件夹下。
    如果对操作系统不是很了解，可以选择方式二

关于程序时间的设置
login.bat 每6个小时（21600 s）自动连一次网（防止学校网络故障）
roboot.bat 一周（7x86400 s, 没办法，windows没有sleep，timeout最大上限只有99999）启动一次（防止连teamviewer久了电脑死机）
如果有人觉得时隔设得太长太短可自行设置

