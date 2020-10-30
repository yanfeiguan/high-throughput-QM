import psycopg2
import time
import socket

from worker import

dbparams = {
    'dbname': 'nmrtest',
    'user': 'nmruser',
    'password': '123456',
    'host': 'localhost',
    'port': 5432,
}


#from bde.gaussian import GaussianRunner


def run_optimization():

    with psycopg2.connect(**dbparams) as conn:
        with conn.cursor() as cur:
            cur.execute("""
            WITH cte AS (
            SELECT id, smiles
            FROM compound
            WHERE status = 'initialized'
            ORDER BY id
            LIMIT 1
            FOR UPDATE
            )

            UPDATE compound SET status = 'running',
                                queued_at = CURRENT_TIMESTAMP,
                                node = %s
            FROM cte
            WHERE compound.id = cte.id
            RETURNING compound.id, compound.smiles;
            """, (socket.gethostname(),))

            cid, smiles = cur.fetchone()
            
    conn.close()

    try:
        
        runner = GaussianRunner(smiles, cid)
        molstr, enthalpy, freeenergy, scfenergy, log = runner.process()
        
        with psycopg2.connect(**dbparams) as conn:
            with conn.cursor() as cur:   
                cur.execute("""
                UPDATE compound
                SET status = 'finished', 
                    mol = %s, enthalpy = %s, 
                    freeenergy = %s, scfenergy= %s,
                    run_at = CURRENT_TIMESTAMP,
                    logfile = %s
                WHERE id = %s;""", 
                (molstr, enthalpy, freeenergy, scfenergy, log, cid))

        conn.close()


    except Exception as ex:
        with psycopg2.connect(**dbparams) as conn:
            with conn.cursor() as cur:
                cur.execute("""
                UPDATE compound
                SET status = 'error', 
                    error = %s, 
                    run_at = CURRENT_TIMESTAMP
                WHERE id = %s;""", (str(ex), cid))  
        
        conn.close()
                
    return cid
                

if __name__ == "__main__":
    
    start_time = time.time()
    run_optimization()

    '''
    # Add a random delay to avoid race conditions at the start of the job
    time.sleep(random.uniform(0, 1*60))
    
    while (time.time() - start_time) < (86400 * 9):  # Time in days
        

        try:
            run_optimization()
            
        except psycopg2.OperationalError:
            time.sleep(5 + random.uniform(0, 60))
    '''