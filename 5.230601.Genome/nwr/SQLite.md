```shell

# 创建数据库，创建一个名为mydatabase.db的SQLite数据库文件
sqlite3 mydatabase.db


# 在SQLite中，数据存储在表格中。下面是一个创建一个名为mytable的表格，包含id、name和age三个列的示例：
CREATE TABLE mytable (
    id INTEGER PRIMARY KEY,
    name TEXT,
    age INTEGER
);

# 插入数据，使用INSERT INTO语句将数据插入表格中：

INSERT INTO mytable (name, age) VALUES ('John', 25);
INSERT INTO mytable (name, age) VALUES ('Jane', 30);

# 查询数据，使用SELECT语句从表格中检索数据：
SELECT * FROM mytable;
SELECT name, age FROM mytable WHERE age > 25;

# 更新数据：












```