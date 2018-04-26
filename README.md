# KRY-RSA
Project for KRY subject
# Usage
`make` compile source code

`make clean` remove binary file

`./kry -g B` generate  `p q n e d` [see here for more information](https://en.wikipedia.org/wiki/RSA_(cryptosystem)) `B` is decimal number defines used bits for representation of `n`

`./kry -e E N M` cipher input message `M` using pair `(E, N)`

`./kry -d D N C` decipher ciphered message `C` using pair `(D, N)`

`./kry -b N` factorize `N` using [Brent-Pollard pho method](http://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf)

P.S. all capital letters(except `B`) are supposed to be hexadecimal numbers, `B` is decimal. 
